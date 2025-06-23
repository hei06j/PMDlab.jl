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
        if !grounded[idx]
            if vmax[idx] < Inf
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 <= vmax[idx]^2)
            end
            if vmin[idx] > 0.0
                JuMP.@constraint(pm.model, vr[t]^2+vi[t]^2 >= vmin[idx]^2)
            end
        end
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
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_negative_sequence(pm::PMD.AbstractUnbalancedACRModel, nw::Int, bus_id::Int, vmnegmax::Real)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmnegsqr)
        # PMD.var(pm, PMD.nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, PMD.nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)

    # real and imaginary components of U-
    vreneg = JuMP.@expression(pm.model,
        (vr_a + a2re*vr_b - a2im*vi_b + are*vr_c - aim*vi_c)/3
    )
    vimneg = JuMP.@expression(pm.model,
        (vi_a + a2re*vi_b + a2im*vr_b + are*vi_c + aim*vr_c)/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@expression(pm.model, vreneg^2+vimneg^2)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmnegsqr <= vmnegmax^2)

    PMD.sol(pm, nw, :bus, bus_id)[:vmnegsqr] = vmnegsqr
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_positive_sequence(pm::PMD.AbstractUnbalancedACRModel, nw::Int, bus_id::Int, vmposmax::Real; vmposmin::Real=0.0)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmpossqr)
        PMD.var(pm, PMD.nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        # PMD.var(pm, PMD.nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@expression(pm.model,
        (vr_a + are*vr_b - aim*vi_b + a2re*vr_c - a2im*vi_c)/3
    )
    vimpos = JuMP.@expression(pm.model,
        (vi_a + are*vi_b + aim*vr_b + a2re*vi_c + a2im*vr_c)/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@expression(pm.model, vrepos^2+vimpos^2)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmpossqr <= vmposmax^2)
    !iszero(vmposmin) ? JuMP.@constraint(pm.model, vmpossqr >= vmposmin^2) : nothing

    PMD.sol(pm, nw, :bus, bus_id)[:vmpossqr] = vmpossqr
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_zero_sequence(pm::PMD.AbstractUnbalancedACRModel, nw::Int, bus_id::Int, vmzeromax::Real)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmzerosqr)
        PMD.var(pm, PMD.nw_id_default)[:vmzerosqr] = Dict{Int, Any}()
    end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    # real and imaginary components of U+
    vrezero = JuMP.@expression(pm.model,
        (vr_a + vr_b + vr_c)/3
    )
    vimzero = JuMP.@expression(pm.model,
        (vi_a + vi_b + vi_c)/3
    )
    # square of magnitude of U+, |U+|^2
    vmzerosqr = JuMP.@expression(pm.model, vrezero^2+vimzero^2)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmzerosqr <= vmzeromax^2)

    PMD.sol(pm, nw, :bus, bus_id)[:vmzerosqr] = vmzerosqr
end



"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_vuf(pm::PMD.AbstractUnbalancedACRModel, nw::Int, bus_id::Int, vufmax::Real)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmpossqr)
        PMD.var(pm, PMD.nw_id_default)[:vmpossqr] = Dict{Int, Any}()
    end
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmnegsqr)
        PMD.var(pm, PMD.nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vuf)
        PMD.var(pm, PMD.nw_id_default)[:vuf] = Dict{Int, Any}()
    end
    (vr_a, vr_b, vr_c) = [PMD.var(pm, nw, :vr, bus_id)[i] for i in 1:3]
    (vi_a, vi_b, vi_c) = [PMD.var(pm, nw, :vi, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    # real and imaginary components of U+
    vrepos = JuMP.@expression(pm.model,
        (vr_a + are*vr_b - aim*vi_b + a2re*vr_c - a2im*vi_c)/3
    )
    vimpos = JuMP.@expression(pm.model,
        (vi_a + are*vi_b + aim*vr_b + a2re*vi_c + a2im*vr_c)/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@expression(pm.model, vrepos^2+vimpos^2)
    # real and imaginary components of U-
    vreneg = JuMP.@expression(pm.model,
        (vr_a + a2re*vr_b - a2im*vi_b + are*vr_c - aim*vi_c)/3
    )
    vimneg = JuMP.@expression(pm.model,
        (vi_a + a2re*vi_b + a2im*vr_b + are*vi_c + aim*vr_c)/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@expression(pm.model, vreneg^2+vimneg^2)
    vuf = JuMP.@expression(pm.model, vreneg / vimneg)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmnegsqr <= vufmax^2*vmpossqr)
    # DEBUGGING: save references for post check
    #PMD.var(pm, PMD.nw_id_default, :vmpossqr)[bus_id] = vmpossqr
    #PMD.var(pm, PMD.nw_id_default, :vmnegsqr)[bus_id] = vmnegsqr

    PMD.sol(pm, nw, :bus, bus_id)[:vuf] = vuf
end


"""
Bounds the current magnitude at both from and to side of a branch
`cr[f_idx]^2 + ci[f_idx]^2 <= c_rating_a^2`
`cr[t_idx]^2 + ci[t_idx]^2 <= c_rating_a^2`
"""
function constraint_mc_current_limit(pm::PMD.AbstractUnbalancedACRModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, c_rating_a::Vector{<:Real})
    (l, f_bus, t_bus) = f_idx
    t_idx = (l, t_bus, f_bus)

    crf =  [PMD.var(pm, nw, :cr, f_idx)[c] for c in f_connections]
    cif =  [PMD.var(pm, nw, :ci, f_idx)[c] for c in f_connections]

    crt =  [PMD.var(pm, nw, :cr, t_idx)[c] for c in t_connections]
    cit =  [PMD.var(pm, nw, :ci, t_idx)[c] for c in t_connections]

    JuMP.@constraint(pm.model, crf.^2 + cif.^2 .<= c_rating_a.^2)
    JuMP.@constraint(pm.model, crt.^2 + cit.^2 .<= c_rating_a.^2)
end
