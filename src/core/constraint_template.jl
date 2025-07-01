"""
    function constraint_mc_voltage_reference(
        pm::AbstractUnbalancedPowerModel,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
    )

Imposes suitable constraints for the voltage at the reference bus
"""
function constraint_mc_voltage_reference(pm::PMD.AbstractUnbalancedPowerModel, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = PMD.ref(pm, nw, :bus, id)
    voltage_source_type = bus["voltage_source_type"]

    if voltage_source_type == "3vm 3va fix"
        constraint_mc_3vm3va_fixed(pm, nw, id, bus["vm"], bus["va"])
    elseif voltage_source_type == "3va fix"
        constraint_mc_3va_fixed(pm, nw, id, bus["va"])
    elseif voltage_source_type == "va fix va diff"
        constraint_mc_vafix_vadiff(pm, nw, id, bus["va"][1], bus["vadelta"])
    elseif voltage_source_type == "va fix seq"
        constraint_mc_vafix_seq(pm, nw, id, bus["va"][1], bus["vm_seq_pos_min"], bus["vm_seq_pos_max"], bus["vm_seq_neg_max"])
    elseif voltage_source_type == "va fix"
        constraint_mc_vafix(pm, nw, id, bus["va"][1])
    end
end


"""
	function constraint_mc_3vm3va_fixed(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_3vm3va_fixed(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, vm::Vector{<:Real}, va::Vector{<:Real})
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    JuMP.@constraint(pm.model, vr .== vm .* cos.(va))
    JuMP.@constraint(pm.model, vi .== vm .* sin.(va))
end


"""
	function constraint_mc_3vm3va_fixed(
		pm::AbstractUnbalancedACPModel,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_3vm3va_fixed(pm::PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, vm::Vector{<:Real}, va::Vector{<:Real})
    JuMP.@constraint(pm.model, PMD.var(pm, nw, :vm, i) .== vm)
    JuMP.@constraint(pm.model, PMD.var(pm, nw, :va, i) .== va)
end

"""
	function constraint_mc_3va_fixed(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		va::Vector{<:Real},
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_3va_fixed(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, va::Vector{<:Real})
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)
    
    for (t, va_ref) in enumerate(va)
        if va_ref == pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] >= 0)
        elseif va_ref == -pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] <= 0)
        elseif va_ref == 0
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        elseif va_ref == pi
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        else
            JuMP.@constraint(pm.model, vi[t] == tan(va_ref)*vr[t])
            # va_ref also implies a sign for vr, vi
            if 0<=va_ref && va_ref <= pi
                JuMP.@constraint(pm.model, vi[t] >= 0)
            else
                JuMP.@constraint(pm.model, vi[t] <= 0)
            end
        end
    end
end

"""
	function constraint_mc_3va_fixed(
		pm::AbstractUnbalancedACPModel,
		nw::Int,
		i::Int,
		va::Vector{<:Real},
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_3va_fixed(pm::PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, va::Vector{<:Real})
    JuMP.@constraint(pm.model, PMD.var(pm, nw, :va, i) .== va)
end



"""
	function constraint_mc_vafix_vadiff(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		va::Real,
        vadelta::Real
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
This function fixes the phase `a` voltage angle.
"""
function constraint_mc_vafix_vadiff(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, va::Real, vadelta::Real)
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)
    
    Re_vavb = JuMP.@expression(pm.model, - (vr[1]*vr[2] + vi[1]*vi[2])/2 + (-vr[1]*vi[2] + vi[1]*vr[2])*sqrt(3)/2 )
    Re_vbvc = JuMP.@expression(pm.model, - (vr[2]*vr[3] + vi[2]*vi[3])/2 + (-vr[2]*vi[3] + vi[2]*vr[3])*sqrt(3)/2 )
    Re_vcva = JuMP.@expression(pm.model, - (vr[3]*vr[1] + vi[3]*vi[1])/2 + (-vr[3]*vi[1] + vi[3]*vr[1])*sqrt(3)/2 )
    
    Im_vavb = JuMP.@expression(pm.model, - (vr[1]*vr[2] + vi[1]*vi[2])*sqrt(3)/2 - (-vr[1]*vi[2] + vi[1]*vr[2])/2 )
    Im_vbvc = JuMP.@expression(pm.model, - (vr[2]*vr[3] + vi[2]*vi[3])*sqrt(3)/2 - (-vr[2]*vi[3] + vi[2]*vr[3])/2 )
    Im_vcva = JuMP.@expression(pm.model, - (vr[3]*vr[1] + vi[3]*vi[1])*sqrt(3)/2 - (-vr[3]*vi[1] + vi[3]*vr[1])/2 )

    if va == pi/2
        JuMP.@constraint(pm.model, vr[1] == 0)
        JuMP.@constraint(pm.model, vi[1] >= 0)
    elseif va == -pi/2
        JuMP.@constraint(pm.model, vr[1] == 0)
        JuMP.@constraint(pm.model, vi[1] <= 0)
    elseif va == 0
        JuMP.@constraint(pm.model, vr[1] >= 0)
        JuMP.@constraint(pm.model, vi[1] == 0)
    elseif va == pi
        JuMP.@constraint(pm.model, vr[1] >= 0)
        JuMP.@constraint(pm.model, vi[1] == 0)
    else
        JuMP.@constraint(pm.model, vi[1] == tan(va)*vr[1])
        # va also implies a sign for vr, vi
        if 0<=va && va <= pi
            JuMP.@constraint(pm.model, vi[1] >= 0)
        else
            JuMP.@constraint(pm.model, vi[1] <= 0)
        end
    end

    JuMP.@constraint(pm.model, Im_vavb <= tan(vadelta * pi/180) * Re_vavb)
    JuMP.@constraint(pm.model, Im_vbvc <= tan(vadelta * pi/180) * Re_vbvc)
    JuMP.@constraint(pm.model, Im_vcva <= tan(vadelta * pi/180) * Re_vcva)

    JuMP.@constraint(pm.model, Im_vavb >= -tan(vadelta * pi/180) * Re_vavb)
    JuMP.@constraint(pm.model, Im_vbvc >= -tan(vadelta * pi/180) * Re_vbvc)
    JuMP.@constraint(pm.model, Im_vcva >= -tan(vadelta * pi/180) * Re_vcva)
end

"""
	function constraint_mc_vafix_vadiff(
		pm::AbstractUnbalancedACPModel,
		nw::Int,
		i::Int,
		va::Real,
        vadelta::Real
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
This function fixes the phase `a` voltage angle.
"""
function constraint_mc_vafix_vadiff(pm::PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, va::Real, vadelta::Real)
    JuMP.@constraint(pm.model, PMD.var(pm, nw, :va, i)[1] == va)

    JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[1] - PMD.var(pm, nw, :va, i)[2] - 120 * pi/180) <= tan((vadelta) * pi/180))
    JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[1] - PMD.var(pm, nw, :va, i)[2] - 120 * pi/180) >= tan((-vadelta) * pi/180))

    JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[2] - PMD.var(pm, nw, :va, i)[3] - 120 * pi/180) <= tan((vadelta) * pi/180))
    JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[2] - PMD.var(pm, nw, :va, i)[3] - 120 * pi/180) >= tan((-vadelta) * pi/180))

    JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[3] - PMD.var(pm, nw, :va, i)[1] - 120 * pi/180) <= tan((vadelta) * pi/180))
    JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[3] - PMD.var(pm, nw, :va, i)[1] - 120 * pi/180) >= tan((-vadelta) * pi/180))

    ## TODO we may need to add another constraint to ensure angles are valid
    # JuMP.@constraint(pm.model, PMD.var(pm, nw, :va, i)[1] + PMD.var(pm, nw, :va, i)[2] + PMD.var(pm, nw, :va, i)[3] == 0)
end


"""
	function constraint_mc_vafix_seq(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		va::Real,
        vm_seq_pos_min::Real,
        vm_seq_pos_max::Real,
        vm_seq_neg_max::Real
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_vafix_seq(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, va::Real, vm_seq_pos_min::Real, vm_seq_pos_max::Real, vm_seq_neg_max::Real)
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)
    
    if va == pi/2
        JuMP.@constraint(pm.model, vr[1] == 0)
        JuMP.@constraint(pm.model, vi[1] >= 0)
    elseif va == -pi/2
        JuMP.@constraint(pm.model, vr[1] == 0)
        JuMP.@constraint(pm.model, vi[1] <= 0)
    elseif va == 0
        JuMP.@constraint(pm.model, vr[1] >= 0)
        JuMP.@constraint(pm.model, vi[1] == 0)
    elseif va == pi
        JuMP.@constraint(pm.model, vr[1] >= 0)
        JuMP.@constraint(pm.model, vi[1] == 0)
    else
        JuMP.@constraint(pm.model, vi[1] == tan(va)*vr[1])
        # va also implies a sign for vr, vi
        if 0<=va && va <= pi
            JuMP.@constraint(pm.model, vi[1] >= 0)
        else
            JuMP.@constraint(pm.model, vi[1] <= 0)
        end
    end

    constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, i, vm_seq_neg_max)
    constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, i, vm_seq_pos_max, vm_seq_pos_min)
end

"""
	function constraint_mc_vafix_seq(
		pm::AbstractUnbalancedACPModel,
		nw::Int,
		i::Int,
		va::Real,
        vm_seq_pos_min::Real,
        vm_seq_pos_max::Real,
        vm_seq_neg_max::Real
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
This function fixes the phase `a` voltage angle.
"""
function constraint_mc_vafix_seq(pm::PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, va::Real, vm_seq_pos_min::Real, vm_seq_pos_max::Real, vm_seq_neg_max::Real)
    JuMP.@constraint(pm.model, PMD.var(pm, nw, :va, i)[1] == va)

    constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, i, vm_seq_neg_max)
    constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, i, vm_seq_pos_max, vm_seq_pos_min)
end



"""
	function constraint_mc_vafix(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		va::Real
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_vafix(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, va::Real)
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)
    
    if va == pi/2
        JuMP.@constraint(pm.model, vr[1] == 0)
        JuMP.@constraint(pm.model, vi[1] >= 0)
    elseif va == -pi/2
        JuMP.@constraint(pm.model, vr[1] == 0)
        JuMP.@constraint(pm.model, vi[1] <= 0)
    elseif va == 0
        JuMP.@constraint(pm.model, vr[1] >= 0)
        JuMP.@constraint(pm.model, vi[1] == 0)
    elseif va == pi
        JuMP.@constraint(pm.model, vr[1] >= 0)
        JuMP.@constraint(pm.model, vi[1] == 0)
    else
        JuMP.@constraint(pm.model, vi[1] == tan(va)*vr[1])
        # va also implies a sign for vr, vi
        if 0<=va && va <= pi
            JuMP.@constraint(pm.model, vi[1] >= 0)
        else
            JuMP.@constraint(pm.model, vi[1] <= 0)
        end
    end
end

"""
	function constraint_mc_vafix(
		pm::AbstractUnbalancedACPModel,
		nw::Int,
		i::Int,
		va::Real
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
This function fixes the phase `a` voltage angle.
"""
function constraint_mc_vafix(pm::PMD.AbstractUnbalancedACPModel, nw::Int, i::Int, va::Real)
    JuMP.@constraint(pm.model, PMD.var(pm, nw, :va, i)[1] == va)
end





"""
    constraint_mc_bus_voltage_balance(pm::AbstractUnbalancedPowerModel, bus_id::Int; nw=nw_id_default)::Nothing

Template function for bus voltage balance constraints.
"""
function constraint_mc_bus_voltage_balance(pm::PMD.AbstractUnbalancedACPModel, bus_id::Int; nw=PMD.nw_id_default)
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
        # PMD.constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"], bus["vm_seq_pos_min"])
        constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"], bus["vm_seq_pos_min"])
        
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
function constraint_mc_bus_voltage_balance(pm::PMD.AbstractUnbalancedACRModel, bus_id::Int; nw=PMD.nw_id_default)
    # @assert(length(PMD.ref(pm, nw, :conductor_ids))==3)

    bus = PMD.ref(pm, nw, :bus, bus_id)

    if haskey(bus, "vm_vuf_max")
        constraint_mc_bus_voltage_magnitude_vuf(pm, nw, bus_id, bus["vm_vuf_max"])
    end

    if haskey(bus, "vm_seq_neg_max")
        constraint_mc_bus_voltage_magnitude_negative_sequence(pm, nw, bus_id, bus["vm_seq_neg_max"])
    end

    if haskey(bus, "vm_seq_pos_max")
        constraint_mc_bus_voltage_magnitude_positive_sequence(pm, nw, bus_id, bus["vm_seq_pos_max"], bus["vm_seq_pos_min"])
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


function constraint_mc_bus_voltage_angle_difference(pm::PMD.AbstractUnbalancedACPModel, i::Int; nw=PMD.nw_id_default)
    # @assert(length(PMD.ref(pm, nw, :conductor_ids))==3)
    bus = PMD.ref(pm, nw, :bus, i)

    if haskey(bus, "va_delta")
        vadelta = bus["va_delta"]

        JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[1] - PMD.var(pm, nw, :va, i)[2] - 120 * pi/180) <= tan((vadelta) * pi/180))
        JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[1] - PMD.var(pm, nw, :va, i)[2] - 120 * pi/180) >= tan((-vadelta) * pi/180))

        JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[2] - PMD.var(pm, nw, :va, i)[3] - 120 * pi/180) <= tan((vadelta) * pi/180))
        JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[2] - PMD.var(pm, nw, :va, i)[3] - 120 * pi/180) >= tan((-vadelta) * pi/180))

        JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[3] - PMD.var(pm, nw, :va, i)[1] - 120 * pi/180) <= tan((vadelta) * pi/180))
        JuMP.@constraint(pm.model, tan(PMD.var(pm, nw, :va, i)[3] - PMD.var(pm, nw, :va, i)[1] - 120 * pi/180) >= tan((-vadelta) * pi/180))
    end
end

function constraint_mc_bus_voltage_angle_difference(pm::PMD.AbstractUnbalancedACRModel, i::Int; nw=PMD.nw_id_default)
    # @assert(length(PMD.ref(pm, nw, :conductor_ids))==3)

    bus = PMD.ref(pm, nw, :bus, i)

    if haskey(bus, "va_delta")
        vadelta = bus["va_delta"]

        vr = PMD.var(pm, nw, :vr, i)
        vi = PMD.var(pm, nw, :vi, i)
        
        Re_vavb = JuMP.@expression(pm.model, - (vr[1]*vr[2] + vi[1]*vi[2])/2 + (-vr[1]*vi[2] + vi[1]*vr[2])*sqrt(3)/2 )
        Re_vbvc = JuMP.@expression(pm.model, - (vr[2]*vr[3] + vi[2]*vi[3])/2 + (-vr[2]*vi[3] + vi[2]*vr[3])*sqrt(3)/2 )
        Re_vcva = JuMP.@expression(pm.model, - (vr[3]*vr[1] + vi[3]*vi[1])/2 + (-vr[3]*vi[1] + vi[3]*vr[1])*sqrt(3)/2 )
        
        Im_vavb = JuMP.@expression(pm.model, - (vr[1]*vr[2] + vi[1]*vi[2])*sqrt(3)/2 - (-vr[1]*vi[2] + vi[1]*vr[2])/2 )
        Im_vbvc = JuMP.@expression(pm.model, - (vr[2]*vr[3] + vi[2]*vi[3])*sqrt(3)/2 - (-vr[2]*vi[3] + vi[2]*vr[3])/2 )
        Im_vcva = JuMP.@expression(pm.model, - (vr[3]*vr[1] + vi[3]*vi[1])*sqrt(3)/2 - (-vr[3]*vi[1] + vi[3]*vr[1])/2 )

        JuMP.@constraint(pm.model, Im_vavb <= tan(vadelta * pi/180) * Re_vavb)
        JuMP.@constraint(pm.model, Im_vbvc <= tan(vadelta * pi/180) * Re_vbvc)
        JuMP.@constraint(pm.model, Im_vcva <= tan(vadelta * pi/180) * Re_vcva)

        JuMP.@constraint(pm.model, Im_vavb >= -tan(vadelta * pi/180) * Re_vavb)
        JuMP.@constraint(pm.model, Im_vbvc >= -tan(vadelta * pi/180) * Re_vbvc)
        JuMP.@constraint(pm.model, Im_vcva >= -tan(vadelta * pi/180) * Re_vcva)
    end
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