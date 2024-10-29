using Pkg
Pkg.activate("./")
using PMDlab
using PowerModelsDistribution
using Ipopt
using JuMP
using LinearAlgebra
optimizer = Ipopt.Optimizer

const PMD = PowerModelsDistribution


## read and parse network data
file = "data/three-wire/network_1/Feeder_2/Master.dss"

eng3w = parse_file(file, transformations=[transform_loops!])
eng3w["settings"]["sbase_default"] = 1        # Change power base here
eng3w["voltage_source"]["source"]["rs"] *= 0  # remove voltage source internal impedance
eng3w["voltage_source"]["source"]["xs"] *= 0  # remove voltage source internal impedance

math3w = transform_data_model(eng3w, kron_reduce=true, phase_project=false)

### changing ref_bus voltage bounds
ref_bus = [i for (i,bus) in math3w["bus"] if occursin("source", bus["name"])]
math3w["bus"]["$(ref_bus[1])"]["vmin"] *= 0.98
math3w["bus"]["$(ref_bus[1])"]["vmax"] *= 1.02

### chaning some of the input data
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


## run optimal power flow AC polar
result3w_acp = solve_mc_opf(math3w, ACPUPowerModel, optimizer)

## run optimal power flow AC rectangular
result3w_acr = solve_mc_opf(math3w, ACRUPowerModel, optimizer)

## run optimal power flow IV rectangular
add_start_vrvi!(math3w)
result3w_ivr = solve_mc_opf(math3w, IVRUPowerModel, optimizer)
