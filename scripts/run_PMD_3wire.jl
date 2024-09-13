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
eng3w["settings"]["sbase_default"] = 1
math3w = transform_data_model(eng3w, kron_reduce=true, phase_project=false)

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

for (d, load) in math3w["load"]
    gen_counter = 2
    if mod(load["index"], 4) == 1
        math3w["gen"]["$gen_counter"] = deepcopy(math3w["gen"]["1"])
        math3w["gen"]["$gen_counter"]["name"] = "$gen_counter"
        math3w["gen"]["$gen_counter"]["index"] = gen_counter
        math3w["gen"]["$gen_counter"]["cost"] = 0.5*math3w["gen"]["1"]["cost"]
        math3w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
        math3w["gen"]["$gen_counter"]["pmax"] = 4*ones(3)
        math3w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
        math3w["gen"]["$gen_counter"]["connections"] = [1;2;3]
        gen_counter = gen_counter + 1
    end
end


### run optimal power flow
add_start_vrvi!(math3w)
result3w = solve_mc_opf(math3w, IVRUPowerModel, optimizer)
