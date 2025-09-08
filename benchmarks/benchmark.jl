
using Comonicon
using DelimitedFiles
using LinearAlgebra
using PMDlab
using PowerModelsDistribution
using Ipopt
using JuMP

import Base: Iterators

using HSL_jll

const PMD = PowerModelsDistribution
const DATA_DIR = joinpath(@__DIR__, "..", "data", "three-wire")
const RESULTS_DIR = joinpath(@__DIR__, "..", "results")

PowerModelsDistribution.silence!()

const VSOURCE_MODEL = [
    "3vm 3va fix", "3va fix", "va fix", "va fix va diff",  "va fix seq",
]

@kwdef struct OptionsBenchmark
    sbase::Float64 = 1.0
    cost_multiplier::Float64 = 1000.0
    line_current_rating::Bool = true
    reduce_lines::Bool = true
    vsource_model::Int = 4
    source_va_rotation::String = "pos"
    initialize_rotation::String = "pos"
    bus_angle_diff_bounds::Bool = false
    Vsequence_bounds::Bool = false
    balanced_impedance::Bool = false
    formulation::String = "ACP"
end

function flag(opt::OptionsBenchmark)
    f1 = opt.line_current_rating ? 1 : 0
    f2 = opt.reduce_lines ? 1 : 0
    f2p = opt.vsource_model

    f3 = opt.source_va_rotation
    f4 = opt.initialize_rotation
    f5 = opt.bus_angle_diff_bounds ? 1 : 0
    f6 = opt.Vsequence_bounds ? 1 : 0
    f7 = opt.balanced_impedance ? 1 : 0

    return "$(opt.formulation)_$(f1)_$(f2)_$(f2p)_$(f3)_$(f4)_$(f5)_$(f6)_$(f7)"
end

function NLPModels.jac_nln_structure!(model::ExaModels.ExaModel, rows, cols)
    NLPModels.jac_structure!(model, rows, cols)
    return rows, cols
end

#=
    JuMP utils
=#

function _num_constraints(model::JuMP.Model)
    non_nl_constraints = sum(JuMP.num_constraints(model, ft, st) for (ft, st) in JuMP.list_of_constraint_types(model) if ft != JuMP.VariableRef)
    return JuMP.num_nonlinear_constraints(model) + non_nl_constraints
end

function _total_callback_time(model)
    # nlp_block = MOI.get(JuMP.unsafe_backend(model), MOI.NLPBlock())
    moi_backend = JuMP.backend(model)
    nlp_block = moi_backend.optimizer.model.nlp_data
    total_callback_time = if isa(nlp_block, MOI.NLPBlockData) && isa(nlp_block.evaluator, MOI.Nonlinear.Evaluator)
        nlp_block.evaluator.eval_objective_timer +
        nlp_block.evaluator.eval_objective_gradient_timer +
        nlp_block.evaluator.eval_constraint_timer +
        nlp_block.evaluator.eval_constraint_jacobian_timer +
        nlp_block.evaluator.eval_hessian_lagrangian_timer
    else
        0.0
    end
    return total_callback_time
end


#=
    Data
=#

function scan_instances(data_dir)
    instances = Tuple{String, String}[]
    for network in readdir(data_dir)
        for feeder in readdir(joinpath(data_dir, network))
            if startswith(feeder, "Feeder")
                push!(instances, (network, feeder))
            end
        end
    end
    return instances
end

function import_data(instance, options)
    eng3w = parse_file(instance, transformations=[transform_loops!])
    PMDlab.augment_eng_3wire!(eng3w; line_current_rating=options.line_current_rating, reduce_lines=options.reduce_lines, sbase=options.sbase)
    math3w = transform_data_model(eng3w, kron_reduce=true, phase_project=true)
    PMDlab.augment_math_3wire!(
        math3w;
        vsource_model=VSOURCE_MODEL[options.vsource_model],
        source_va_rotation=options.source_va_rotation,
        bus_angle_diff_bounds=options.bus_angle_diff_bounds,
        Vsequence_bounds=options.Vsequence_bounds,
        balanced_impedance=options.balanced_impedance,
        initialize_rotation=options.initialize_rotation,
        cost_multiplier=options.cost_multiplier,
    )  # changing some of the input data
    return math3w
end

function build_pm_model(instance, options)
    type = if options.formulation == "ACP"
        ACPUPowerModel
    elseif options.formulation == "IVR"
        IVRUPowerModel
    elseif options.formulation == "ACR"
        ACRUPowerModel
    end
    return PMD.instantiate_mc_model(
        import_data(instance, options),
        type,
        PMD.build_mc_opf;
    )
end

#=
    Solvers
=#
function solve_ipopt(model)
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_attribute(model, "max_iter", 1000)
    JuMP.set_attribute(model, "max_wall_time", 600.0)
    JuMP.set_attribute(model, "hsllib", HSL_jll.libhsl_path)
    JuMP.set_attribute(model, "linear_solver", "ma27")
    solve_time = @elapsed begin
        JuMP.optimize!(model)
    end
    success = JuMP.termination_status(model) == MOI.LOCALLY_SOLVED
    return (
        nvar=JuMP.num_variables(model),
        ncon=_num_constraints(model),
        success=success,
        objective=JuMP.objective_value(model),
        solve_time=solve_time,
        callback_time=_total_callback_time(model),
        iter=MOI.get(model, MOI.BarrierIterations()),
    )
end

#=
    Benchmark
=#

function run_benchmark(solver, instances, options)
    nexp = length(instances)
    results = zeros(nexp, 7)

    k = 1
    for (network, feeder) in instances
        @info "$(network) $(feeder)"
        pm = build_pm_model(joinpath(DATA_DIR, network, feeder, "Master.dss"), options)
        res = solver(pm.model)

        results[k, 1] = res.nvar
        results[k, 2] = res.ncon
        results[k, 3] = res.success
        results[k, 4] = res.objective
        results[k, 5] = res.iter
        results[k, 6] = res.solve_time
        results[k, 7] = res.callback_time
        k += 1
    end

    c1, c2 = [r[1] for r in instances], [r[2] for r in instances]
    return [c1 c2 results]
end

function run_all_benchmarks(solver, instances, dump_dir)
    if !isdir(joinpath(RESULTS_DIR, dump_dir))
        mkdir(joinpath(RESULTS_DIR, dump_dir))
    end

    itr1 = ["ACP", "IVR", "ACR"]
    itr2 = 1:5
    itr3 = ["pos", "neg", "zero"]
    itr4 = [true, false]
    itr5 = [true, false]
    itr6 = [true, false]

    iterators = Iterators.product(itr1, itr2, itr3, itr4, itr5, itr6)
    for (form, vsource_model, svarot, seqbnd, bangle, bimp) in iterators
        try
            options = OptionsBenchmark(
                formulation=form,
                vsource_model=vsource_model,
                source_va_rotation=svarot,
                initialize_rotation=svarot,
                Vsequence_bounds=seqbnd,
                bus_angle_diff_bounds=bangle,
                balanced_impedance=bimp,
            )
            results = run_benchmark(solver, instances, options)
            dump_file = joinpath(RESULTS_DIR, dump_dir, "$(flag(options)).txt")
            writedlm(dump_file, results)
        catch ex
            println(ex)
            continue
        end
    end
end

function main()
    instances = scan_instances(DATA_DIR)
    run_all_benchmarks(solve_ipopt, instances, "ipopt")
end

main()

