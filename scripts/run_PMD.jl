using Pkg
Pkg.activate("./")
using PMDlab

using PowerModelsDistribution
using Ipopt
using JuMP
using LinearAlgebra
optimizer = Ipopt.Optimizer

const PMD = PowerModelsDistribution


##
# file = "./data/opendss/TestFeeder/1/Master.dss"
file = "data/enwl_4w/network_1/Feeder_2/Master.dss"
eng4w = parse_file(file, transformations=[transform_loops!,remove_all_bounds!])
eng4w["settings"]["sbase_default"] = 1
math4w = transform_data_model(eng4w, kron_reduce=false, phase_project=false)
add_start_vrvi!(math4w)
# fix voltage source 

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

vsource_branch, vsource_bus = PMDlab.find_voltage_source_branch_bus(math4w)
transformer_primary_bus = PMDlab.find_transformer_primary_bus(math4w)
math4w["gen"]["1"]["gen_bus"] = transformer_primary_bus
math4w["bus"]["$transformer_primary_bus"] = deepcopy(math4w["bus"]["$vsource_bus"])
math4w["bus"]["$transformer_primary_bus"]["bus_i"] = transformer_primary_bus
math4w["bus"]["$transformer_primary_bus"]["index"] = transformer_primary_bus
delete!(math4w["branch"], vsource_branch)
delete!(math4w["bus"], "$vsource_bus")


pm4w = instantiate_mc_model(
        math4w,
        IVRENPowerModel,
        build_mc_opf_alt_4w;
        # ref_extensions=ref_extensions,
        multinetwork=false,
        # kwargs...
    )


result4w = solve_mc_opf(math4w, IVRENPowerModel, optimizer)
for (i, bus) in result4w["solution"]["bus"]
    # bus["vm"] = hypot.(bus["vr"], bus["vi"])
    if length(bus["vr"]) >= 4
        bus["vmn"] = abs.((bus["vr"][1:3] .+ im.*bus["vi"][1:3]) .- (bus["vr"][4] + im.*bus["vi"][4]))
    end
end



## phase to neutral transform 

# file3wpn = "./data/opendss/TestFeeder/2/Master.dss"
file3wpn = pathof(PMDlab)[1:end-22]*"data/opendss/ENWL/enwl_phase_to_neutral_embedded/network_1/Feeder_2/Master.dss"
eng3wpn = parse_file(file3wpn, transformations=[remove_all_bounds!])
eng3wpn["is_kron_reduced"] = true
eng3wpn["settings"]["sbase_default"] = 1

# math = transform_data_model(eng, kron_reduce=false, phase_project=false)
math3wpn = transform_data_model(eng3wpn, kron_reduce=false, phase_project=false)
vsource_branch, vsource_bus = find_voltage_source_branch_bus(math3wpn)
math3wpn["branch"][vsource_branch]
math3wpn["conductor_ids"] = math3wpn["conductor_ids"][1:3]


for (i,bus) in math3wpn["bus"]
    if haskey(bus, "grounded")
        bus["grounded"] = bus["grounded"][1:3]
    end
    if haskey(bus, "terminals")
        bus["terminals"] = bus["terminals"][1:3]
    end 
    if haskey(bus, "vmax")
        bus["vmax"] = bus["vmax"][1:3]
    end
    if haskey(bus, "vmin")
        bus["vmin"] = bus["vmin"][1:3]
    end
    bus["vmin"] = [0.9;0.9;0.9]
    bus["vmax"] = [1.1;1.1;1.1]
end

for (l,branch) in math3wpn["branch"]
    if haskey(branch, "t_connections")
        branch["t_connections"] = branch["t_connections"][1:3]
    end
    if haskey(branch, "f_connections")
        branch["f_connections"] = branch["f_connections"][1:3]
    end
    if haskey(branch, "br_r")
        branch["br_r"] = branch["br_r"][1:3,1:3]
    end
    if haskey(branch, "br_x")
        branch["br_x"] = branch["br_x"][1:3,1:3]
    end
    if haskey(branch, "g_to")
        branch["g_to"] = branch["g_to"][1:3,1:3]
    end
    if haskey(branch, "g_fr")
        branch["g_fr"] = branch["g_fr"][1:3,1:3]
    end
    if haskey(branch, "b_to")
        branch["b_to"] = branch["b_to"][1:3,1:3]
    end
    if haskey(branch, "b_fr")
        branch["b_fr"] = branch["b_fr"][1:3,1:3]
    end
    if haskey(branch, "c_rating_a")
        branch["c_rating_a"] = branch["c_rating_a"][1:3]
    end
end

for (t,transformer) in math3wpn["transformer"]
    if haskey(transformer, "t_connections")
        transformer["t_connections"] = transformer["t_connections"][1:3]
    end
    if haskey(transformer, "f_connections")
        transformer["f_connections"] = transformer["f_connections"][1:3]
    end
end

transformer_primary_bus = find_transformer_primary_bus(math3wpn)
math3wpn["gen"]["1"]["gen_bus"] = transformer_primary_bus
math3wpn["bus"]["$transformer_primary_bus"] = deepcopy(math3wpn["bus"]["$vsource_bus"])
math3wpn["bus"]["$transformer_primary_bus"]["bus_i"] = transformer_primary_bus
math3wpn["bus"]["$transformer_primary_bus"]["index"] = transformer_primary_bus
delete!(math3wpn["branch"], vsource_branch)
delete!(math3wpn["bus"], "$vsource_bus")

math3wpn["bus"]["$transformer_primary_bus"]["vmin"] = [1;1;1.0]
math3wpn["bus"]["$transformer_primary_bus"]["vmax"] = [1;1;1.0]



for (g,gen) in math3wpn["gen"]
    gen["connections"] = gen["connections"][1:3]
    gen["vg"] = gen["vg"][1:3]
    gen["pg"] = gen["pg"][1:3]
    gen["qg"] = gen["qg"][1:3]
    gen["pmax"] = gen["pmax"][1:3]
    gen["pmin"] = gen["pmin"][1:3]
    gen["qmax"] = gen["qmax"][1:3]
    gen["qmin"] = gen["qmin"][1:3]
    gen["cost"] = 1000 .* gen["cost"]
end

for (d, load) in math3wpn["load"]
    gen_counter = 2
    if mod(load["index"], 4) == 1
        math3wpn["gen"]["$gen_counter"] = deepcopy(math3wpn["gen"]["1"])
        math3wpn["gen"]["$gen_counter"]["name"] = "$gen_counter"
        math3wpn["gen"]["$gen_counter"]["index"] = gen_counter
        math3wpn["gen"]["$gen_counter"]["cost"] = 0.5*math3wpn["gen"]["1"]["cost"]
        math3wpn["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
        math3wpn["gen"]["$gen_counter"]["pmax"] = 4*ones(3)
        math3wpn["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
        math3wpn["gen"]["$gen_counter"]["connections"] = [1;2;3]
        gen_counter = gen_counter + 1
    end
end

for (l,load) in math3wpn["load"]
    load["connections"] = load["connections"][1:1]
end

include("../scripts/3winitialization.jl")
add_start_voltage_3w!(math3wpn, coordinates=:rectangular)
# add_start_vrvi!(math3wpn)

pm3wpn = instantiate_mc_model(
    math3wpn,
        IVRUPowerModel ,
        PMDlab.build_mc_opf_alt_3w;
        # ref_extensions=ref_extensions,
        multinetwork=false,
        # kwargs...
    )

# print(pm3wpn.model)

result3wpn = solve_mc_opf(math3wpn, IVRUPowerModel, optimizer)
for (i, bus) in result3wpn["solution"]["bus"]
    bus["vmn"] = hypot.(bus["vr"], bus["vi"])
end



##
@show (result4w["solve_time"])
@show (result3wpn["solve_time"])
@show (result4w["objective"])
@show (result3wpn["objective"])

@show (result4w["solution"]["bus"]["2"])
@show (result3wpn["solution"]["bus"]["2"])

@show (result4w["solution"]["bus"]["621"])
@show (result3wpn["solution"]["bus"]["621"])

@show (result4w["solution"]["gen"]["1"]["pg"])
@show (result3wpn["solution"]["gen"]["1"]["pg"])

@show (result4w["solution"]["gen"]["2"]["pg"])
@show (result3wpn["solution"]["gen"]["2"]["pg"])


# @show math4w["branch"]["1"]["br_r"]
# @show math3wpn["branch"]["1"]["br_r"]



##
