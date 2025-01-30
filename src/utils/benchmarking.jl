function augment_eng_3wire!(eng; line_current_rating=true, sbase=1)

    if line_current_rating
        add_linecode_normaps!(eng)           # add branch current limit
    end

    eng["settings"]["sbase_default"] = sbase    # change power base here
    eng["voltage_source"]["source"]["rs"] *= 0  # remove voltage source internal impedance
    eng["voltage_source"]["source"]["xs"] *= 0  # remove voltage source internal impedance

end


function augment_math_3wire!(math; relax_vsource_vm=true, Vsequence_bounds=true, cost_multiplier=1000)
    math["is_kron_reduced"] = true

    fix_data_parsing_issues!(math)       # ensuring length-3 vector across 3-wire data model

    if Vsequence_bounds
        add_sequence_voltage_bounds!(math)   # add sequence voltage bounds
    end


    ### changing ref_bus voltage bounds
    if relax_vsource_vm
        ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
        math["bus"]["$(ref_bus[1])"]["vmin"] *= 0.98
        math["bus"]["$(ref_bus[1])"]["vmax"] *= 1.02
    end

    for (i,bus) in math["bus"]
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            # bus["grounded"] .=  0
        end
    end
    
    for (g,gen) in math["gen"]
        gen["cost"] = gen["cost"] * cost_multiplier
    end
    
    gen_counter = length(math["gen"])
    for (d, load) in math["load"]
        if mod(load["index"], 4) == 1
            gen_counter += 1
            math["gen"]["$gen_counter"] = deepcopy(math["gen"]["1"])
            math["gen"]["$gen_counter"]["name"] = "$gen_counter"
            math["gen"]["$gen_counter"]["index"] = gen_counter
            math["gen"]["$gen_counter"]["cost"] = 0.5*math["gen"]["1"]["cost"]
            math["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
            math["gen"]["$gen_counter"]["pmax"] = 4*ones(3)
            math["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
            math["gen"]["$gen_counter"]["connections"] = [1;2;3]
        end
        ### change every 10th load to constant impedance
        if mod(load["index"], 10) == 1
            load["model"] = PMD.IMPEDANCE
        end
    end
end


function augment_network_data_4wire!(math; relax_vsource_vm=true, cost_multiplier=1000)
    ### changing ref_bus voltage bounds
    if relax_vsource_vm
        ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
        math["bus"]["$(ref_bus[1])"]["vmin"] *= 0.98
        math["bus"]["$(ref_bus[1])"]["vmax"] *= 1.02
    end

    for (i,bus) in math["bus"]
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            bus["grounded"] .=  0
        end
    end

    for (g,gen) in math["gen"]
        gen["cost"] = gen["cost"] * cost_multiplier
    end

    gen_counter = length(math["gen"])
    shunt_counter = length(math["shunt"])
    for (d, load) in math["load"]
        if mod(load["index"], 4) == 1
            gen_counter += 1
            math["gen"]["$gen_counter"] = deepcopy(math["gen"]["1"])
            math["gen"]["$gen_counter"]["name"] = "gen_$gen_counter"
            math["gen"]["$gen_counter"]["index"] = gen_counter
            math["gen"]["$gen_counter"]["cost"] = 0.5*math["gen"]["1"]["cost"]
            math["gen"]["$gen_counter"]["gen_bus"] = load["load_bus"]
            math["gen"]["$gen_counter"]["pmax"] = 4*ones(3)
            math["gen"]["$gen_counter"]["pmin"] = 0.0*ones(3)
            math["gen"]["$gen_counter"]["connections"] = [1;2;3;4]
        end
        ### change every 10th load to constant impedance
        if mod(load["index"], 10) == 1
            load["model"] = PMD.IMPEDANCE
        end
    end
end

function find_voltage_source_branch_bus(math)
    for (b, branch) in math["branch"]
        if branch["source_id"] == "voltage_source.source"
            return b, branch["f_bus"]
        end
    end
    return error()
end


function add_linecode_normaps!(eng)
    linecode_normamps = Dict(
        "lc1" => 75, 
        "lc2" => 107,
        "lc3" => 129,
        "lc4" => 167,
        "lc5" => 223,
        "lc6" => 58,
        "lc7" => 281,
        "lc8" => 419,
        "lc9" => 197)
    for (l, linecode) in eng["linecode"]
        linecode["cm_ub"] = linecode_normamps[l] * ones(length(linecode["cm_ub"]))
    end
end


function add_sequence_voltage_bounds!(math)
    tr_bus_list = []
    if !isempty(math["transformer"])
        tr_bus_list = get_transformer_buses(math)
    end
    
    for (i, bus) in math["bus"]
        if parse(Int, i) ∉ tr_bus_list
            # bus["vm_seq_pos_max"] = 1.1
            bus["vm_seq_neg_max"] = 0.02
            # bus["vm_seq_zero_max"] = 0.02
            # bus["vm_vuf_max"] = 2# 0.02
        end
    end
end


function get_transformer_buses(math)
    bus_list = []
    for (i, tr) in math["transformer"]
        tr_fbus = tr["f_bus"]
        tr_tbus = tr["t_bus"]
        append!(bus_list, tr_fbus)
        append!(bus_list, tr_tbus)

        tr_branches = [b for (b, branch) in math["branch"] if (branch["f_bus"] ∈ [tr_fbus, tr_tbus] || branch["t_bus"] ∈ [tr_fbus, tr_tbus])]
        for b in tr_branches
            append!(bus_list, math["branch"]["$b"]["f_bus"])
            append!(bus_list, math["branch"]["$b"]["t_bus"])
        end
    end
    return unique(bus_list)
end


function fix_data_parsing_issues!(math)
    @assert length(math["conductor_ids"]) == 3 "The mathematical model is not for a three wire case"

    for (i, tr) in math["transformer"]
        if length(tr["f_connections"]) == 4
            tr["f_connections"] = tr["f_connections"][1:3]
        end
        if length(tr["t_connections"]) == 4
            tr["t_connections"] = tr["t_connections"][1:3]
        end
    end

    for (i, gen) in math["gen"]
        if length(gen["connections"]) == 4
            gen["connections"] = gen["connections"][1:3]
        end
    end

    for (i, bus) in math["bus"]
        if occursin("source", bus["name"])
            if length(bus["terminals"]) == 4
                bus["terminals"] = bus["terminals"][1:3]
                bus["grounded"] = bus["grounded"][1:3]
                bus["vmin"] = bus["vmin"][1:3]
                bus["vmax"] = bus["vmax"][1:3]
                bus["vm"] = bus["vm"][1:3]
                bus["va"] = bus["va"][1:3]
            end

        end
    end

    if length(math["conductor_ids"]) == 4
        math["conductor_ids"] = math["conductor_ids"][1:3]
    end
    
end