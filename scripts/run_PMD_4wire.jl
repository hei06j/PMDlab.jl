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
file = "data/four-wire/network_1/Feeder_2/Master.dss"

eng4w = parse_file(file, transformations=[transform_loops!])
eng4w["settings"]["sbase_default"] = 1        # Change power base here
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
        math4w["gen"]["$gen_counter"]["name"] = "gen_$gen_counter"
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

## run optimal power flow  IV rectangular
result4w_ivr = solve_mc_opf(math4w, IVRENPowerModel, optimizer)

for (i, bus) in result4w_ivr["solution"]["bus"]
    if length(bus["vr"]) >= 4
        bus["vmn"] = abs.((bus["vr"][1:3] .+ im.*bus["vi"][1:3]) .- (bus["vr"][4] + im.*bus["vi"][4]))
    end
end

