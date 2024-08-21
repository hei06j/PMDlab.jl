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
file = "data/enwl_4w/network_1/Feeder_2/Master.dss"

eng4w = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
eng4w["settings"]["sbase_default"] = 1
math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)

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

for (d, load) in math4w["load"]
    gen_counter = 2
    if mod(load["index"], 4) == 1
        math4w["gen"]["$gen_counter"] = deepcopy(math4w["gen"]["1"])
        math4w["gen"]["$gen_counter"]["name"] = "$gen_counter"
        math4w["gen"]["$gen_counter"]["index"] = gen_counter
        math4w["gen"]["$gen_counter"]["cost"] = 0.5*math4w["gen"]["1"]["cost"]
        math4w["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
        math4w["gen"]["$gen_counter"]["pmax"] = 4*ones(3)
        math4w["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
        math4w["gen"]["$gen_counter"]["connections"] = [1;2;3;4]
        gen_counter = gen_counter + 1
    end
end


### run optimal power flow
add_start_vrvi!(math4w)
result4w = solve_mc_opf(math4w, IVRENPowerModel, optimizer)


for (i, bus) in result4w["solution"]["bus"]
    if length(bus["vr"]) >= 4
        bus["vmn"] = abs.((bus["vr"][1:3] .+ im.*bus["vi"][1:3]) .- (bus["vr"][4] + im.*bus["vi"][4]))
    end
end
