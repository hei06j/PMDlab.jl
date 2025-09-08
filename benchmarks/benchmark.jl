
using Comonicon
using DelimitedFiles
using PMDlab
using PowerModelsDistribution
using Ipopt
using KNITRO
using NLPModelsKnitro
using JuMP
using LinearAlgebra

using MadNLP
using MadNLPHSL
using ExaModels
using NLPModels

using HSL_jll

struct ThreeWire end
struct FourWire end

const PMD = PowerModelsDistribution
const DATA_DIR = joinpath(@__DIR__, "data")

PowerModelsDistribution.silence!()

function NLPModels.jac_nln_structure!(model::ExaModels.ExaModel, rows, cols)
    NLPModels.jac_structure!(model, rows, cols)
    return rows, cols
end
function NLPModels.jac_nln_coord!(model::ExaModels.ExaModel, x, jac)
    NLPModels.jac_coord!(model, x, jac)
    return jac
end
function NLPModels.hess_structure(model::ExaModels.ExaModel)
    nnzh = NLPModels.get_nnzh(model)
    rows = zeros(Int, nnzh)
    cols = zeros(Int, nnzh)
    NLPModels.hess_structure!(model, rows, cols)
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
            push!(instances, (network, feeder))
        end
    end
    return instances
end

function import_data(instance, ::ThreeWire)
    eng3w = parse_file(instance, transformations=[transform_loops!])
    eng3w["settings"]["sbase"] = 1
    eng3w["voltage_source"]["source"]["rs"] *= 0  # remove voltage source internal impedance
    eng3w["voltage_source"]["source"]["xs"] *= 0  # remove voltage source internal impedance

    math3w = transform_data_model(eng3w, kron_reduce=true, phase_project=false)

    ref_bus = [i for (i,bus) in math3w["bus"] if occursin("source", bus["name"])]
    math3w["bus"]["$(ref_bus[1])"]["vmin"] *= 0.98
    math3w["bus"]["$(ref_bus[1])"]["vmax"] *= 1.02

    ### changing some of the input data
    for (i,bus) in math3w["bus"]
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            # bus["grounded"] .=  0
        end
    end

    for (g,gen) in math3w["gen"]
        gen["cost"] = 1000 .* gen["cost"]
    end

    gen_counter = length(math3w["gen"])
    for (d, load) in math3w["load"]
        if mod(load["index"], 4) == 1
            gen_counter = gen_counter + 1
            math3w["gen"]["$gen_counter"] = deepcopy(math3w["gen"]["1"])
            math3w["gen"]["$gen_counter"]["name"] = "$gen_counter"
            math3w["gen"]["$gen_counter"]["index"] = gen_counter
            math3w["gen"]["$gen_counter"]["cost"] = 0.5*math3w["gen"]["1"]["cost"]
            math3w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
            math3w["gen"]["$gen_counter"]["pmax"] = 4*ones(3)
            math3w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
            math3w["gen"]["$gen_counter"]["connections"] = [1;2;3]
        end
        ### change every 10th load to constant impedance
        if mod(load["index"], 10) == 1
            load["model"] = IMPEDANCE
        end
    end
    data_math = transform_data_model(
        math3w;
    )
    return data_math
end

function import_data(instance, ::FourWire)
    eng4w = parse_file(instance, transformations=[transform_loops!])
    eng4w["settings"]["sbase"] = 1
    eng4w["voltage_source"]["source"]["rs"] *= 0  # remove voltage source internal impedance
    eng4w["voltage_source"]["source"]["xs"] *= 0  # remove voltage source internal impedance

    ### adding grounding to every second load, in the engineering model
    if !haskey(eng4w, "shunt")
        eng4w["shunt"] = Dict{String,Any}()
    end
    shunt_counter = length(eng4w["shunt"])
    for (d, load) in enumerate(eng4w["load"])
        if mod(d, 2) == 1
            shunt_counter += 1
            load_bus = last(load)["bus"]
            eng4w["shunt"]["$shunt_counter"] = Dict{String,Any}()
            eng4w["shunt"]["$shunt_counter"]["source_id"] = "reactor.grounding_load_$d"
            eng4w["shunt"]["$shunt_counter"]["status"] = ENABLED
            eng4w["shunt"]["$shunt_counter"]["connections"] = [4, 5]
            eng4w["shunt"]["$shunt_counter"]["bus"] = load_bus
            eng4w["shunt"]["$shunt_counter"]["gs"] = [0.1  -0.1 ; -0.1   0.1]
            eng4w["shunt"]["$shunt_counter"]["bs"] = [0.0 0.0 ; 0.0 0.0]
            eng4w["bus"]["$load_bus"]["terminals"] = [1,2,3,4,5]
            eng4w["bus"]["$load_bus"]["grounded"] = [5]
            eng4w["bus"]["$load_bus"]["rg"] = [0.0]
            eng4w["bus"]["$load_bus"]["xg"] = [0.0]
        end
    end

    math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)

    ### changing ref_bus voltage bounds
    ref_bus = [i for (i,bus) in math4w["bus"] if occursin("source", bus["name"])]
    math4w["bus"]["$(ref_bus[1])"]["vmin"] *= 0.98
    math4w["bus"]["$(ref_bus[1])"]["vmax"] *= 1.02

    ### chaning some of the input data
    for (i,bus) in math4w["bus"]
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            bus["grounded"] .=  0
        end
    end

    for (g,gen) in math4w["gen"]
        gen["cost"] = 1000 .* gen["cost"]
    end

    gen_counter = length(math4w["gen"])
    shunt_counter = length(math4w["shunt"])
    for (d, load) in math4w["load"]
        if mod(load["index"], 4) == 1
            gen_counter = gen_counter + 1
            math4w["gen"]["$gen_counter"] = deepcopy(math4w["gen"]["1"])
            math4w["gen"]["$gen_counter"]["name"] = "$gen_counter"
            math4w["gen"]["$gen_counter"]["index"] = gen_counter
            math4w["gen"]["$gen_counter"]["cost"] = 0.5*math4w["gen"]["1"]["cost"]
            math4w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
            math4w["gen"]["$gen_counter"]["pmax"] = 4*ones(3)
            math4w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
            math4w["gen"]["$gen_counter"]["connections"] = [1;2;3;4]
        end
        ### change every 10th load to constant impedance
        if mod(load["index"], 10) == 1
            load["model"] = IMPEDANCE
        end
    end
    add_start_vrvi!(math4w)
    data_math = transform_data_model(
        math4w;
    )
    return data_math
end

function build_pm_model(instance, type, wire)
    return instantiate_mc_model(
        import_data(instance, wire),
        type,
        build_mc_opf;
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

function solve_madnlp(model)
    JuMP.set_optimizer(model, MadNLP.Optimizer)
    JuMP.set_attribute(model, "max_iter", 1000)
    JuMP.set_attribute(model, "max_wall_time", 600.0)
    JuMP.set_attribute(model, "linear_solver", Ma27Solver)
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

function solve_knitro(model)
    JuMP.set_optimizer(model, KNITRO.Optimizer)
    JuMP.set_attribute(model, "maxit", 1000)
    JuMP.set_attribute(model, "maxtime_real", 600.0)
    JuMP.set_attribute(model, "linsolver", 4) # ma27
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

function solve_exa_madnlp(model)
    nlp = ExaModels.ExaModel(model)
    total_time = @elapsed begin
        res = madnlp(
            nlp;
            linear_solver=Ma27Solver,
            max_iter=1000,
            max_wall_time=600.0,
        )
    end
    success = res.status == MadNLP.SOLVE_SUCCEEDED
    return (
        nvar=NLPModels.get_nvar(nlp),
        ncon=NLPModels.get_ncon(nlp),
        success=success,
        objective=res.objective,
        solve_time=total_time,
        callback_time=res.counters.eval_function_time,
        iter=res.iter,
    )
end

function solve_exa_knitro(model)
    nlp = ExaModels.ExaModel(model)
    total_time = @elapsed begin
        res = knitro(
            nlp;
            linsolver=4,
            maxit=1000,
            maxtime_real=600.0,
        )
    end
    success = res.status == :first_order
    return (
        nvar=NLPModels.get_nvar(nlp),
        ncon=NLPModels.get_ncon(nlp),
        success=success,
        objective=res.objective,
        solve_time=total_time,
        callback_time=0.0,
        iter=res.iter,
    )
end

#=
    Benchmark
=#

function run_benchmark(solver, instances, formulation, wire)
    nexp = length(instances)
    results = zeros(nexp, 7)
    subdir = isa(wire, ThreeWire) ? "three-wire" : "four-wire"

    k = 1
    for (network, feeder) in instances
        @info "$(network) $(feeder)"
        pm = build_pm_model(joinpath(DATA_DIR, subdir, network, feeder, "Master.dss"), formulation, wire)
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

@main function main(; solver="ipopt", form="acp", nwire="3-wire")
    if nwire == "3-wire"
        instances = scan_instances(joinpath(DATA_DIR, "three-wire"))
        formulation = if form == "acp"
            ACPUPowerModel
        elseif form == "acr"
            ACRUPowerModel
        elseif form == "ivr"
            IVRUPowerModel
        end
        wire = ThreeWire()
        nw = 3
    elseif nwire == "4-wire"
        instances = scan_instances(joinpath(DATA_DIR, "four-wire"))
        form = "ivr"
        formulation = IVRENPowerModel
        wire = FourWire()
        nw = 4
    end

    if solver == "ipopt"
        results = run_benchmark(solve_ipopt, instances, formulation, wire)
        dump_file = joinpath(@__DIR__, "results", "ipopt_ma27_$(nw)w_$(form).txt")
    elseif solver == "knitro"
        results = run_benchmark(solve_knitro, instances, formulation, wire)
        dump_file = joinpath(@__DIR__, "results", "knitro_ma27_$(nw)w_$(form).txt")
    elseif solver == "exaknitro"
        results = run_benchmark(solve_exa_knitro, instances, formulation, wire)
        dump_file = joinpath(@__DIR__, "results", "exaknitro_ma27_$(nw)w_$(form).txt")
    elseif solver == "madnlp"
        results = run_benchmark(solve_madnlp, instances, formulation, wire)
        dump_file = joinpath(@__DIR__, "results", "madnlp_ma27_$(nw)w_$(form).txt")
    elseif solver == "examadnlp"
        results = run_benchmark(solve_exa_madnlp, instances, formulation, wire)
        dump_file = joinpath(@__DIR__, "results", "examadnlp_ma27_$(nw)w_$(form).txt")
    end
    writedlm(dump_file, results)
end

