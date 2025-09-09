using Pkg
Pkg.activate("./")
using PMDlab
using PowerModelsDistribution
using Ipopt
using JuMP
using LinearAlgebra
using DataFrames
using CSV

import JuMP._CONSTRAINT_LIMIT_FOR_PRINTING
JuMP._CONSTRAINT_LIMIT_FOR_PRINTING[] = 10000

optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes", "max_iter"=>3000)
optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>3000)

const PMD = PowerModelsDistribution

PMD.silence!()

## run optimal power flow IV rectangular
# file = "data/three-wire/network_23/Feeder_3/Master.dss"     # three-wire without transformer  SMALLEST network
# file = "data/three-wire/network_1/Feeder_1/Master.dss"      # three-wire without transformer  LARGER network

file = "data/three-wire/network_16/Feeder_2/Master.dss"     # three-wire without transformer  SMALLEST network
file = "data/three-wire/network_8/Feeder_2/Master.dss"     # three-wire without transformer  SMALLEST network
# file = "data/three-wire/network_19/Feeder_1/Master.dss"     # three-wire without transformer  SMALLEST network

""" line_current_rating adds thermal rating to lines.
    OPTIONS: true, false """
line_current_rating = true

""" reduce_lines.
    OPTIONS: true, false """
reduce_lines = true

""" vsource_model is the model of the voltage source bus.
    OPTIONS: "3vm 3va fix", "3va fix", "va fix", "va fix va diff",  "va fix seq" """
# vsource_model = "va fix va diff"
vsource_model = "3vm 3va fix"

""" source_va_rotation initializes the bus voltages with angles in the specified sequence roration for the reference bus.
    OPTIONS: "pos", "neg", "zero" """
source_va_rotation = "pos"

""" initialize_rotation initializes the bus voltages with angles in the specified sequence roration for all buses other than reference bus.
    OPTIONS: "pos", "neg", "zero" """
initialize_rotation = "pos"

""" bus_angle_diff_bounds adds bus voltage angle difference bounds on buses other than the reference bus.
    OPTIONS: true, false """
bus_angle_diff_bounds = false

""" Vsequence_bounds adds voltage sequence bounds on buses other than the reference bus.
    OPTIONS: true, false """
Vsequence_bounds = false

""" balanced_impedance makes all impedances balanced.
    OPTIONS: true, false """
balanced_impedance = false

""" Formulations for OPF
    OPTIONS: "IVR", "ACR", "ACP" """
formulation = "IVR"

"""
if either  bus_angle_diff_bounds or Vsequence_bounds are true, source_va_rotation and initialize_rotation should be the same.
    otherwise, they can be different, but we need to test feasibility and alignment.

if vsource_model is "va fix", "va fix va diff", "va fix seq", source_va_rotation does not matter as it only takes into account the first va,
    so initialize_rotation can be any option, but still check for convergence.
"""


"""
What to test:
- "reduce line" and "line_current_rating" always true
- "vsource_model" test all 5 to check convergence
- "source_va_rotation" and "initialize_rotation" test compatible options (pos-pos, neg-neg, zero-zero)
- "Vsequence_bounds" both cases
- "bus_angle_diff_bounds" both cases
- "balanced_impedance" both cases
- "formulation" all cases - starting with ACP, IVR to test if ACP improves with angle difference bounds
- optionally test "sbase" (constraint scaling) and "cost_multiplier" (objective scaling)
"""
##
eng3w = parse_file(file, transformations=[transform_loops!])
PMDlab.augment_eng_3wire!(eng3w; line_current_rating=line_current_rating, reduce_lines=reduce_lines, sbase=100)
math3w = transform_data_model(eng3w, kron_reduce=true, phase_project=false)
PMDlab.augment_math_3wire!(math3w; vsource_model=vsource_model, source_va_rotation=source_va_rotation, bus_angle_diff_bounds=bus_angle_diff_bounds, Vsequence_bounds=Vsequence_bounds, balanced_impedance=balanced_impedance, initialize_rotation=initialize_rotation, cost_multiplier=1000)  # changing some of the input data


result3w = []
if formulation == "IVR"
    pm = PMD.instantiate_mc_model(math3w, IVRUPowerModel, PMDlab.build_mc_opf);
    result3w = PMD.solve_mc_model(math3w, IVRUPowerModel, optimizer, PMDlab.build_mc_opf)
    vm = [sqrt.(bus["vr"].^2 .+ bus["vi"].^2) for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]
    va = [angle.(bus["vr"] .+ im*bus["vi"]).*180/pi for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]

elseif formulation == "ACR"
    pm = PMD.instantiate_mc_model(math3w, ACRUPowerModel, PMDlab.build_mc_opf);  # "va fix va diff" does not work
    result3w = PMD.solve_mc_model(math3w, ACRUPowerModel, optimizer, PMDlab.build_mc_opf)
    vm = [sqrt.(bus["vr"].^2 .+ bus["vi"].^2) for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]
    va = [angle.(bus["vr"] .+ im*bus["vi"]).*180/pi for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]

elseif formulation == "ACP"
    pm = PMD.instantiate_mc_model(math3w, ACPUPowerModel, PMDlab.build_mc_opf);
    result3w = PMD.solve_mc_model(math3w, ACPUPowerModel, optimizer, PMDlab.build_mc_opf)
    vm = [bus["vm"] for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]
    va = [bus["va"].*180/pi for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]

end



## Some basic checks

load_id_bus = [(i, load["load_bus"]) for (i,load) in math3w["load"]]
load_s_gap = []
for (i,busid) in load_id_bus
    load_connections = math3w["load"]["$i"]["connections"]

    load_c = result3w["solution"]["load"]["$i"]["crd_bus"] .+ im * result3w["solution"]["load"]["$i"]["cid_bus"]
    load_v = result3w["solution"]["bus"]["$busid"]["vr"][load_connections] .+ im * result3w["solution"]["bus"]["$busid"]["vi"][load_connections]
    load_s0 = result3w["solution"]["load"]["$i"]["pd"] .+ im * result3w["solution"]["load"]["$i"]["qd"]
    load_s = load_v .* conj(load_c)
    push!(load_s_gap, load_s .- load_s0)
end
maxmax = maximum(max.([maximum(abs.(vec)) for vec in load_s_gap]))


gen_id_bus = [(i, gen["gen_bus"]) for (i,gen) in math3w["gen"]]
gen_s_gap = []
for (i,busid) in gen_id_bus
    gen_connections = math3w["gen"]["$i"]["connections"]

    gen_c = result3w["solution"]["gen"]["$i"]["crg_bus"] .+ im * result3w["solution"]["gen"]["$i"]["cig_bus"]
    gen_v = result3w["solution"]["bus"]["$busid"]["vr"][gen_connections] .+ im * result3w["solution"]["bus"]["$busid"]["vi"][gen_connections]
    gen_s0 = result3w["solution"]["gen"]["$i"]["pg"] .+ im * result3w["solution"]["gen"]["$i"]["qg"]
    gen_s = gen_v .* conj(gen_c)
    push!(gen_s_gap, gen_s .- gen_s0)
end
maxmax = maximum(max.([maximum(abs.(vec)) for vec in gen_s_gap]))


total_gen = sum(sum(gen["pg"]+im*gen["qg"] for (i,gen) in result3w["solution"]["gen"]))
total_load = sum(sum(load["pd"]+im*load["qd"] for (i,load) in result3w["solution"]["load"]))

total_losses = total_gen .- total_load
loss_percent = total_losses / total_gen * 100

##

vneqseq = [sqrt(bus["vmnegsqr"]) for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]
vposseq = [sqrt(bus["vmpossqr"]) for (i, bus) in result3w["solution"]["bus"] if occursin("source", math3w["bus"]["$i"]["name"])]

pg = [gen["pg"] for (i,gen) in result3w["solution"]["gen"]]



# output_file = "model.json"
# open(output_file, "w") do f
#     # for i in eachindex(PMD.(pm, 0, ))
#     #     println(f, )
#     # end
#     print(f, pm.model)
# end

##

## For all the test cases, the following new constraints are checked:
### - negative sequence bound should be binding in some of the testcase
### - current limits also should be binding in some of the testcase
ENWL_3W_EMBD_DIR = "data/three-wire"

checks_df = DataFrames.DataFrame()

for network_dir in readdir(ENWL_3W_EMBD_DIR)
    mn = match(r"network_(\d+)$"i, network_dir)
    if !isnothing(mn)
        n = mn.captures[1]
        for feeder_dir in readdir("$ENWL_3W_EMBD_DIR/$network_dir")
            mf = match(r"feeder_(\d+)$"i, feeder_dir)
            if !isnothing(mf)
                f = mf.captures[1]
                case_dir = "$ENWL_3W_EMBD_DIR/$network_dir/$feeder_dir"

                eng3w = parse_file(case_dir*"/Master.dss", transformations=[transform_loops!])
                PMDlab.augment_eng_3wire!(eng3w; line_current_rating=line_current_rating, reduce_lines=reduce_lines, sbase=1)
                math3w = transform_data_model(eng3w, kron_reduce=true, phase_project=true)
                PMDlab.augment_math_3wire!(math3w; vsource_model=vsource_model, source_va_rotation=source_va_rotation, bus_angle_diff_bounds=bus_angle_diff_bounds, Vsequence_bounds=Vsequence_bounds, balanced_impedance=balanced_impedance, initialize_rotation=initialize_rotation, cost_multiplier=1000)  # changing some of the input data
                pm = PMD.instantiate_mc_model(math3w, IVRUPowerModel, PMDlab.build_mc_opf);
                result3w = PMD.solve_mc_model(math3w, IVRUPowerModel, optimizer, PMDlab.build_mc_opf)
                df, _ = PMDlab.check_active_bounds(result3w, math3w)

                append!(checks_df, df)
            end
        end
    end
end

CSV.write("IVR_checks.csv", checks_df)