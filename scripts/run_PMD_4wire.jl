using Pkg
Pkg.activate("./")
using PMDlab
using PowerModelsDistribution
using Ipopt
using JuMP
using LinearAlgebra
optimizer = Ipopt.Optimizer

const PMD = PowerModelsDistribution

# function test()

## read and parse network data
file = "data/four-wire/network_1/Feeder_2/Master.dss"
# file = "data/four-wire/network_10/Feeder_5/Master.dss"

eng4w = parse_file(file, transformations=[transform_loops!])
PMDlab.add_linecode_normaps!(eng4w)           # add branch current limit
eng4w["settings"]["sbase_default"] = 1        # change power base here
eng4w["voltage_source"]["source"]["rs"] *= 0  # remove voltage source internal impedance
eng4w["voltage_source"]["source"]["xs"] *= 0  # remove voltage source internal impedance
reduce_line_series!(eng4w, remove_original_lines=true)

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
PMDlab.add_sequence_voltage_bounds!(math4w)     # add sequence voltage bounds


### chaning some of the input data
PMDlab.augment_network_data_4wire!(math4w)

##
add_start_vrvi!(math4w)

## run optimal power flow  IV rectangular
result4w_ivr = solve_mc_opf(math4w, IVRENPowerModel, optimizer)

for (i, bus) in result4w_ivr["solution"]["bus"]
    if length(bus["vr"]) >= 4
        bus["vmn"] = abs.((bus["vr"][1:3] .+ im.*bus["vi"][1:3]) .- (bus["vr"][4] + im.*bus["vi"][4]))
    end
end

