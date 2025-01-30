"""
    constraint_mc_bus_voltage_balance(pm::AbstractUnbalancedPowerModel, bus_id::Int; nw=nw_id_default)::Nothing

Template function for bus voltage balance constraints.
"""
function constraint_mc_bus_voltage_balance(pm::PMD.AbstractUnbalancedACPModel, bus_id::Int; nw=PMD.nw_id_default)::Nothing
    # @assert(length(PMD.ref(pm, nw, :conductor_ids))==3)

    bus = PMD.ref(pm, nw, :bus, bus_id)

    if haskey(bus, "vm_vuf_max")
        # PMD.constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
        constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    end

    if haskey(bus, "vm_seq_neg_max")
        # PMD.constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, bus["vm_seq_neg_max"])
        constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, bus["vm_seq_neg_max"])
    end

    if haskey(bus, "vm_seq_pos_max")
        # PMD.constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"])
        constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"])
    end

    if haskey(bus, "vm_seq_zero_max")
        # PMD.constraint_mc_bus_voltage_magnitude_zero_sequence(pm, nw, bus_id, bus["vm_seq_zero_max"])
        constraint_mc_bus_voltage_magnitude_zero_sequence(pm, nw, bus_id, bus["vm_seq_zero_max"])
    end

    if haskey(bus, "vm_ll_min")|| haskey(bus, "vm_ll_max")
        vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : fill(0, 3)
        vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : fill(Inf, 3)
        PMD.constraint_mc_bus_voltage_magnitude_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    end
    nothing
end


"""
    constraint_mc_bus_voltage_balance(pm::AbstractUnbalancedACRModel, bus_id::Int; nw=nw_id_default)::Nothing

Template function for bus voltage balance constraints.
"""
function constraint_mc_bus_voltage_balance(pm::PMD.AbstractUnbalancedACRModel, bus_id::Int; nw=PMD.nw_id_default)::Nothing
    # @assert(length(PMD.ref(pm, nw, :conductor_ids))==3)

    bus = PMD.ref(pm, nw, :bus, bus_id)

    if haskey(bus, "vm_vuf_max")
        constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    end

    if haskey(bus, "vm_seq_neg_max")
        constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, bus["vm_seq_neg_max"])
    end

    if haskey(bus, "vm_seq_pos_max")
        constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"])
    end

    if haskey(bus, "vm_seq_zero_max")
        constraint_mc_bus_voltage_magnitude_zero_sequence(pm, nw, bus_id, bus["vm_seq_zero_max"])
    end

    if haskey(bus, "vm_ll_min")|| haskey(bus, "vm_ll_max")
        vm_ll_min = haskey(bus, "vm_ll_min") ? bus["vm_ll_min"] : fill(0, 3)
        vm_ll_max = haskey(bus, "vm_ll_max") ? bus["vm_ll_max"] : fill(Inf, 3)
        PMD.constraint_mc_bus_voltage_magnitude_ll(pm, nw, bus_id, vm_ll_min, vm_ll_max)
    end
    nothing
end



# """
# Bounds the current magnitude at both from and to side of a branch
# `cr[f_idx]^2 + ci[f_idx]^2 <= c_rating_a^2`
# `cr[t_idx]^2 + ci[t_idx]^2 <= c_rating_a^2`
# """
# function constraint_mc_current_limit(pm::PMD.AbstractUnbalancedACPModel, i::Int; nw=PMD.nw_id_default)::Nothing
#     branch = _PMD.ref(pm, nw, :branch, i)
#     f_bus = branch["f_bus"]
#     t_bus = branch["t_bus"]
#     f_idx = (i, f_bus, t_bus)

#     constraint_mc_current_limit(pm, nw, f_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"])
# end


"""
Bounds the current magnitude at both from and to side of a branch
`cr[f_idx]^2 + ci[f_idx]^2 <= c_rating_a^2`
`cr[t_idx]^2 + ci[t_idx]^2 <= c_rating_a^2`
"""
function constraint_mc_current_limit(pm::PMD.AbstractUnbalancedACRModel, i::Int; nw=PMD.nw_id_default)
    branch = PMD.ref(pm, nw, :branch, i)
    f_bus = branch["f_bus"]
    t_bus = branch["t_bus"]
    f_idx = (i, f_bus, t_bus)
    
    if all(branch["c_rating_a"] .< Inf)
        constraint_mc_current_limit(pm, nw, f_idx, branch["f_connections"], branch["t_connections"], branch["c_rating_a"])
    end
end