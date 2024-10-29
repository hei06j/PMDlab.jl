
function run_mc_opf(math, formulation, post, optimizer)
    pm = PMD.instantiate_mc_model(
        math,
        formulation,
        post;
        # ref_extensions=ref_extensions,
        multinetwork=false,
        # kwargs...
    )
    result = PMD.solve_mc_opf(math, formulation, optimizer)
    return result
end


"""
    function constraint_mc_voltage_reference(
        pm::ExplicitNeutralModels,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
    )

Imposes suitable constraints for the voltage at the reference bus
"""
function constraint_mc_voltage_reference(pm::PMD.AbstractUnbalancedPowerModel, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = PMD.ref(pm, nw, :bus, id)
    terminals = bus["terminals"]
    grounded = bus["grounded"]

    if haskey(bus, "va") && !haskey(bus, "vm")
        constraint_mc_theta_ref(pm, id; nw=nw)
    elseif haskey(bus, "vm") && !haskey(bus, "va")
        constraint_mc_voltage_magnitude_fixed(pm, nw, id, bus["vm"], terminals, grounded)
    elseif haskey(bus, "vm") && haskey(bus, "va")
        constraint_mc_voltage_fixed(pm, nw, id, bus["vm"], bus["va"], terminals, grounded)
    end
end


"""
    constraint_mc_theta_ref(pm::AbstractUnbalancedPowerModel, i::Int; nw::Int=nw_id_default)::Nothing

Template function for reference angle constraints.
"""
function constraint_mc_theta_ref(pm::PMD.AbstractUnbalancedPowerModel, i::Int; nw::Int=PMD.nw_id_default)::Nothing
    bus = PMD.ref(pm, nw, :bus, i)
    terminals = bus["terminals"]
    if haskey(bus, "va")
        va_ref = get(PMD.ref(pm, nw, :bus, i), "va", [deg2rad.([0.0, -120.0, 120.0])..., zeros(length(terminals))...][terminals])
        constraint_mc_theta_ref(pm, nw, i, va_ref)
    end
    nothing
end


"Creates phase angle constraints at reference buses"
function constraint_mc_theta_ref(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, va_ref::Vector{<:Real})
    terminals = PMD.ref(pm, nw, :bus, i, "terminals")
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    # deal with cases first where tan(theta)==Inf or tan(theta)==0
    for (idx, t) in enumerate(terminals)
        if va_ref[t] == pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] >= 0)
        elseif va_ref[t] == -pi/2
            JuMP.@constraint(pm.model, vr[t] == 0)
            JuMP.@constraint(pm.model, vi[t] <= 0)
        elseif va_ref[t] == 0
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        elseif va_ref[t] == pi
            JuMP.@constraint(pm.model, vr[t] >= 0)
            JuMP.@constraint(pm.model, vi[t] == 0)
        else
            JuMP.@constraint(pm.model, vi[t] == tan(va_ref[idx])*vr[t])
            # va_ref also implies a sign for vr, vi
            if 0<=va_ref[t] && va_ref[t] <= pi
                JuMP.@constraint(pm.model, vi[t] >= 0)
            else
                JuMP.@constraint(pm.model, vi[t] <= 0)
            end
        end
    end
end



"""
	function constraint_mc_voltage_magnitude_fixed(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_voltage_magnitude_fixed(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, vm::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    @assert iszero(vm[grounded]) "Infeasible model; the voltage magnitude of a grounded terminal is fixed to a non-zero value."
    idxs = findall((!).(grounded))

    JuMP.@constraint(pm.model, [i in idxs], vr[terminals[i]].^2 + vi[terminals[i]].^2 == vm[i].^2)
end


"""
	function constraint_mc_voltage_fixed(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		vm::Vector{<:Real},
		va::Vector{<:Real},
		terminals::Vector{Int},
		grounded::Vector{Bool}
	)

Fixes the voltage variables at bus `i` to `vm.*exp.(im*va)`
"""
function constraint_mc_voltage_fixed(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, vm::Vector{<:Real}, va::Vector{<:Real}, terminals::Vector{Int}, grounded::Vector{Bool})
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    idxs = findall((!).(grounded))

    JuMP.@constraint(pm.model, [i in idxs], vr[terminals[i]]==vm[i]*cos(va[i]))
    JuMP.@constraint(pm.model, [i in idxs], vi[terminals[i]]==vm[i]*sin(va[i]))
end


"""
    function constraint_mc_voltage_absolute(
        pm::AbstractUnbalancedACRModel,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
    )

Imposes absolute voltage magnitude bounds for models with explicit neutrals
"""
function constraint_mc_voltage_absolute(pm::PMD.AbstractUnbalancedACRModel, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = PMD.ref(pm, nw, :bus, id)

    constraint_mc_voltage_absolute(pm, nw, id, bus["terminals"], bus["grounded"], bus["vmin"], bus["vmax"])
end


"""
    function constraint_mc_voltage_pairwise(
        pm::AbstractUnbalancedACRModel,
        id::Int;
        nw::Int=nw_id_default,
        bounded::Bool=true,
        report::Bool=true,
    )

Imposes pairwise voltage magnitude bounds, i.e. magnitude bounds on the voltage between to terminals, for models with explicit neutrals
"""
function constraint_mc_voltage_pairwise(pm::PMD.AbstractUnbalancedACRModel, id::Int; nw::Int=PMD.nw_id_default, bounded::Bool=true, report::Bool=true)
    bus = PMD.ref(pm, nw, :bus, id)

    vm_pair_lb = bus["vm_pair_lb"]
    vm_pair_ub = bus["vm_pair_ub"]

    constraint_mc_voltage_pairwise(pm, nw, id, vm_pair_lb, vm_pair_ub)
end


"""
	function constraint_mc_voltage_absolute(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		vmin::Vector{<:Real},
		vmax::Vector{<:Real};
		report::Bool=true
	)

Imposes absolute voltage magnitude bounds for models with explicit neutrals
"""
function constraint_mc_voltage_absolute(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, terminals::Vector{Int}, grounded::Vector{Bool}, vmin::Vector{<:Real}, vmax::Vector{<:Real}; report::Bool=true)
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    # ungrounded_terminals = terminals[(!).(grounded)]
    for (idx,t) in enumerate(terminals)
        # if !grounded[idx]
            if vmax[idx] < Inf
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 <= vmax[idx]^2)
            end
            if vmin[idx] > 0.0
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 >= vmin[idx]^2)
            end
        # end
    end
end



"""
	function constraint_mc_voltage_pairwise(
		pm::AbstractUnbalancedACRModel,
		nw::Int,
		i::Int,
		terminals::Vector{Int},
		grounded::Vector{Bool},
		vm_pair_lb::Vector,
		vm_pair_ub::Vector;
		report::Bool=true
	)

Imposes pairwise voltage magnitude bounds, i.e. magnitude bounds on the voltage between to terminals, for models with explicit neutrals
"""
function constraint_mc_voltage_pairwise(pm::PMD.AbstractUnbalancedACRModel, nw::Int, i::Int, vm_pair_lb::Vector{<:Tuple{Any,Any,Real}}, vm_pair_ub::Vector{<:Tuple{Any,Any,Real}}; report::Bool=true)
    vr = PMD.var(pm, nw, :vr, i)
    vi = PMD.var(pm, nw, :vi, i)

    for (a,b,lb) in vm_pair_lb
        if lb > 0.0
            JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 >= lb^2)
        end
    end

    for (a,b,ub) in vm_pair_ub
        if ub < Inf
            JuMP.@constraint(pm.model, (vr[a]-vr[b])^2 + (vi[a]-vi[b])^2 <= ub^2)
        end
    end
end


"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""

"""
	function build_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

constructor for OPF in current-voltage variable space
"""
function build_mc_opf(pm::PMD.AbstractUnbalancedIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_switch_current(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_load_current(pm)

    # Constraints
    # for i in PMD.ids(pm, :ref_buses)
    #     PMD.constraint_mc_theta_ref(pm, i)
    # end
    # Constraints
    for i in ids(pm, :bus)

        if i in ids(pm, :ref_buses)
            constraint_mc_voltage_reference(pm, i)
        end

        constraint_mc_voltage_absolute(pm, i)
        constraint_mc_voltage_pairwise(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.onstraint_mc_current_balance(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)

        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_current_limit(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i)
    end

    # Objective
    PMD.objective_mc_min_fuel_cost(pm)
end


# function build_mc_opf_alt_4w(pm::PMD.AbstractExplicitNeutralIVRModel)
#     # Variables
#     PMD.variable_mc_bus_voltage(pm)
#     PMD.variable_mc_branch_current(pm)
#     PMD.variable_mc_load_current(pm)
#     PMD.variable_mc_load_power(pm)
#     PMD.variable_mc_generator_current(pm)
#     PMD.variable_mc_generator_power(pm)
#     PMD.variable_mc_transformer_current(pm)
#     PMD.variable_mc_transformer_power(pm)
#     PMD.variable_mc_switch_current(pm)

#     # Constraints
#     for i in PMD.ids(pm, :bus)

#         if i in PMD.ids(pm, :ref_buses)
#             PMD.constraint_mc_voltage_reference(pm, i)
#         end

#         PMD.constraint_mc_voltage_absolute(pm, i)
#         PMD.constraint_mc_voltage_pairwise(pm, i)
#     end

#     # components should be constrained before KCL, or the bus current variables might be undefined

#     for id in PMD.ids(pm, :gen)
#         PMD.constraint_mc_generator_power(pm, id)
#         PMD.constraint_mc_generator_current(pm, id)
#     end

#     for id in PMD.ids(pm, :load)
#         PMD.constraint_mc_load_power(pm, id)
#         PMD.constraint_mc_load_current(pm, id)
#     end

#     for i in PMD.ids(pm, :transformer)
#         PMD.constraint_mc_transformer_voltage(pm, i)
#         PMD.constraint_mc_transformer_current(pm, i)

#         PMD.constraint_mc_transformer_thermal_limit(pm, i)
#     end

#     for i in PMD.ids(pm, :branch)
#         PMD.constraint_mc_current_from(pm, i)
#         PMD.constraint_mc_current_to(pm, i)
#         PMD.constraint_mc_bus_voltage_drop(pm, i)

#         PMD.constraint_mc_branch_current_limit(pm, i)
#         # PMD.constraint_mc_thermal_limit_from(pm, i)
#         # PMD.constraint_mc_thermal_limit_to(pm, i)
#     end

#     # for i in PMD.ids(pm, :switch)
#     #     PMD.constraint_mc_switch_current(pm, i)
#     #     PMD.constraint_mc_switch_state(pm, i)

#     #     PMD.constraint_mc_switch_current_limit(pm, i)
#     #     PMD.constraint_mc_switch_thermal_limit(pm, i)
#     # end

#     for i in PMD.ids(pm, :bus)
#         PMD.constraint_mc_current_balance(pm, i)
#     end

#     # Objective
#     PMD.objective_mc_min_fuel_cost(pm)
# end


