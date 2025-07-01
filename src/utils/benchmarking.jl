function make_impedances_symmetric!(math)
    for (l, branch) in math["branch"]
        r = branch["br_r"]
        x = branch["br_x"]
        branch["br_r"] .= (r .+ r[[2, 3, 1],[2, 3, 1]] + r[[3, 1, 2],[3, 1,2]])./3
        branch["br_x"] .= (x .+ x[[2, 3, 1],[2, 3, 1]] + x[[3, 1, 2],[3, 1,2]])./3
    end
end

function augment_eng_3wire!(eng; line_current_rating=true, reduce_lines=true, sbase=1)

    if line_current_rating
        add_linecode_normaps!(eng)           # add branch current limit
    end

    if reduce_lines
        PMD.reduce_line_series!(eng, remove_original_lines=true)
    end

    eng["settings"]["sbase_default"] = sbase    # change power base here
    eng["voltage_source"]["source"]["rs"] *= 0  # remove voltage source internal impedance
    eng["voltage_source"]["source"]["xs"] *= 0  # remove voltage source internal impedance

end


function augment_math_3wire!(math; vsource_model="3vm 3va fix", source_va_rotation="pos", bus_angle_diff_bounds=true, Vsequence_bounds=true, balanced_impedance=false, initialize_rotation="pos", cost_multiplier=1000)
    math["is_kron_reduced"] = true
    
    fix_data_parsing_issues!(math)       # ensuring length-3 vector across 3-wire data model

    if balanced_impedance
        make_impedances_symmetric!(math)     # make impedance matrices balanced
    end

    if bus_angle_diff_bounds || Vsequence_bounds
        @assert source_va_rotation == initialize_rotation "source va rotation and initialization are inconsistent in presence of bounds"
    end

    if source_va_rotation !== initialize_rotation
        @warn("source va rotation and initialization are inconsistent")
    end

    # Vsequence_bounds = "pos", "neg", "zero", "vuf" TODO
    if Vsequence_bounds 
        @assert source_va_rotation == "pos"  "The bus voltage sequence bounds should only be active for positive sequence initialization"
        add_sequence_voltage_bounds!(math)   # add sequence voltage bounds
    end

    if bus_angle_diff_bounds 
        @assert source_va_rotation == "pos"  "The bus voltage angle difference bounds should only be active for positive sequence initialization"
        add_bus_angle_diff_bounds!(math)   # add sequence voltage bounds
    end

    ### source_va_rotation = "pos", "neg", "zero"
    va_vector = []
    if source_va_rotation == "pos"
        va_vector = [0, -120, 120] .* pi/180
    elseif source_va_rotation == "neg"
        va_vector = [0, 120, -120] .* pi/180
    elseif source_va_rotation == "zero"
        va_vector = [0, 0, 0] .* pi/180
    else
        Error("source_va_rotation option not identified")
    end

    ### vsource_model = "3vm 3va fix", "3va fix", "va fix va diff", "va fix seq"
    if vsource_model == "3va fix"
        ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
        math["bus"]["$(ref_bus[1])"]["vmin"] = 0.9 .* ones(3)
        math["bus"]["$(ref_bus[1])"]["vmax"] = 1.1 .* ones(3)
        # math["bus"]["$(ref_bus[1])"]["vm"] = [1, 1, 1]
        math["bus"]["$(ref_bus[1])"]["va"] = va_vector
        math["bus"]["$(ref_bus[1])"]["voltage_source_type"] = "3va fix"
    elseif vsource_model == "va fix va diff"
        ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
        math["bus"]["$(ref_bus[1])"]["vmin"] .= 0   # 0.9 .* ones(3)
        math["bus"]["$(ref_bus[1])"]["vmax"] .= Inf # 1.1 .* ones(3)
        # math["bus"]["$(ref_bus[1])"]["vm"] = [1, 1, 1]
        math["bus"]["$(ref_bus[1])"]["va"] = va_vector
        math["bus"]["$(ref_bus[1])"]["vadelta"] = 30
        math["bus"]["$(ref_bus[1])"]["voltage_source_type"] = "va fix va diff"
    elseif vsource_model == "va fix"
        ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
        math["bus"]["$(ref_bus[1])"]["vmin"] .= 0   # 0.9 .* ones(3)
        math["bus"]["$(ref_bus[1])"]["vmax"] .= Inf # 1.1 .* ones(3)
        # math["bus"]["$(ref_bus[1])"]["vm"] = [1, 1, 1]
        math["bus"]["$(ref_bus[1])"]["va"] = va_vector
        math["bus"]["$(ref_bus[1])"]["voltage_source_type"] = "va fix"
    elseif vsource_model == "va fix seq"
        ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
        math["bus"]["$(ref_bus[1])"]["vmin"] .= 0   # 0.9 .* ones(3)
        math["bus"]["$(ref_bus[1])"]["vmax"] .= Inf # 1.1 .* ones(3)
        math["bus"]["$(ref_bus[1])"]["vm_seq_pos_min"] = 0.9
        math["bus"]["$(ref_bus[1])"]["vm_seq_pos_max"] = 1.1
        math["bus"]["$(ref_bus[1])"]["vm_seq_neg_max"] = 0.02
        # math["bus"]["$(ref_bus[1])"]["vm_seq_zero_max"] = 0.02
        # math["bus"]["$(ref_bus[1])"]["vm_vuf_max"] = 2# 0.02
        # math["bus"]["$(ref_bus[1])"]["vm"] = [1, 1, 1]
        math["bus"]["$(ref_bus[1])"]["va"] = va_vector
        math["bus"]["$(ref_bus[1])"]["voltage_source_type"] = "va fix seq"
    else # "3vm 3va fix"
        ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
        math["bus"]["$(ref_bus[1])"]["vmin"] .= 0
        math["bus"]["$(ref_bus[1])"]["vmax"] .= Inf
        math["bus"]["$(ref_bus[1])"]["vm"] = [1, 1, 1]
        math["bus"]["$(ref_bus[1])"]["va"] = va_vector
        math["bus"]["$(ref_bus[1])"]["voltage_source_type"] = "3vm 3va fix"
    end
    # if vsource_model_vm
    #     ref_bus = [i for (i,bus) in math["bus"] if occursin("source", bus["name"])]
    #     math["bus"]["$(ref_bus[1])"]["vmin"] *= 0.98
    #     math["bus"]["$(ref_bus[1])"]["vmax"] *= 1.02
    # end

    for (i,bus) in math["bus"]
        if bus["bus_type"] != 3 && !startswith(bus["source_id"], "transformer")
            # bus["vm_pair_lb"] = [(1, 4, 0.9);(2, 4, 0.9);(3, 4, 0.9)]
            # bus["vm_pair_ub"] = [(1, 4, 1.1);(2, 4, 1.1);(3, 4, 1.1)]
            # bus["grounded"] .=  0
            bus["vmin"] = 0.9 .* [1, 1, 1]
            bus["vmax"] = 1.1 .* [1, 1, 1]
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


    ## how to initialise
    if initialize_rotation == "pos"
        add_start_positive_sequence!(math)  # this works with bus angle difference and voltage sequence constraints, and source_va_rotation="pos"

    elseif initialize_rotation == "neg"
        add_start_negative_sequence!(math)  # this works without bus angle difference and voltage sequence constraints, and source_va_rotation="neg"
        
    elseif initialize_rotation == "zero"
        add_start_zero_sequence!(math)      # this works without bus angle difference and voltage sequence constraints, and source_va_rotation="zero"
    else
        PMD.add_start_voltage!(math, coordinates=:rectangular, explicit_neutral=false)
    end
end


function augment_network_data_4wire!(math; vsource_model_vm=true, cost_multiplier=1000)
    ### changing ref_bus voltage bounds
    if vsource_model_vm
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
            # bus["vm_seq_pos_min"] = 0.8
            bus["vm_seq_neg_max"] = 0.02
            # bus["vm_seq_zero_max"] = 0.02
            # bus["vm_vuf_max"] = 2# 0.02
        end
    end
end

function add_bus_angle_diff_bounds!(math)
    tr_bus_list = []
    if !isempty(math["transformer"])
        tr_bus_list = get_transformer_buses(math)
    end
    
    for (i, bus) in math["bus"]
        if parse(Int, i) ∉ tr_bus_list
            bus["va_delta"] = 30  # degrees
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


### initialization functions
function add_start_positive_sequence!(math)
    vrvi_start = [1.0 + 0.0im, -0.4999999999999998 - 0.8660254037844387im, -0.4999999999999998 + 0.8660254037844387im]
    for (i, bus) in math["bus"]
        bus["vr_start"] = real.(vrvi_start)
        bus["vr_start"] = imag.(vrvi_start)
        bus["vm_start"] = abs.(vrvi_start)
        bus["va_start"] = angle.(vrvi_start)
    end
end

function add_start_negative_sequence!(math)
    vrvi_start = [1.0 + 0.0im, -0.4999999999999998 + 0.8660254037844387im, -0.4999999999999998 - 0.8660254037844387im]
    for (i, bus) in math["bus"]
        bus["vr_start"] = real.(vrvi_start)
        bus["vr_start"] = imag.(vrvi_start)
        bus["vm_start"] = abs.(vrvi_start)
        bus["va_start"] = angle.(vrvi_start)
    end
end

function add_start_zero_sequence!(math)
    vrvi_start = [1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im]
    for (i, bus) in math["bus"]
        bus["vr_start"] = real.(vrvi_start)
        bus["vr_start"] = imag.(vrvi_start)
        bus["vm_start"] = abs.(vrvi_start)
        bus["va_start"] = angle.(vrvi_start)
    end
end