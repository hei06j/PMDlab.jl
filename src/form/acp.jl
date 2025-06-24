"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_negative_sequence(pm::PMD.AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vmnegmax::Real)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmnegsqr)
        PMD.var(pm, PMD.nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMD.var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [PMD.var(pm, nw, :va, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)

    ### real and imaginary components of U-
    vreneg = JuMP.NonlinearExpr[]
    vimneg = JuMP.NonlinearExpr[]
    vmnegsqr = JuMP.NonlinearExpr[]
    vreneg = JuMP.@expression(pm.model,
        (vm_a*cos(va_a) + a2re*vm_b*cos(va_b) - a2im*vm_b*sin(va_b) + are*vm_c*cos(va_c) - aim*vm_c*sin(va_c))/3
    )
    vimneg = JuMP.@expression(pm.model,
        (vm_a*sin(va_a) + a2re*vm_b*sin(va_b) + a2im*vm_b*cos(va_b) + are*vm_c*sin(va_c) + aim*vm_c*cos(va_c))/3
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
function constraint_mc_bus_voltage_magnitude_positive_sequence(pm::PMD.AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vmposmax::Real, vmposmin::Real)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmpossqr)
        PMD.var(pm, PMD.nw_id_default)[:vmpossqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMD.var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [PMD.var(pm, nw, :va, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)
    ### real and imaginary components of U+
    # vrepos = JuMP.NonlinearExpr[]
    # vimpos = JuMP.NonlinearExpr[]
    # vmpossqr = JuMP.NonlinearExpr[]
    # vrepos = JuMP.@expression(pm.model,
    #     (vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3
    # )
    # vimpos = JuMP.@expression(pm.model,
    #     (vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3
    # )

    # square of magnitude of U+, |U+|^2
    # vmpossqr = JuMP.@expression(pm.model, vrepos^2+vimpos^2)
    # finally, apply constraint
    # JuMP.@NLconstraint(pm.model, vrepos^2+vimpos^2 <= vmposmax^2)
    # JuMP.@NLconstraint(pm.model, vrepos^2+vimpos^2 >= vmposmin^2)

    JuMP.@constraint(pm.model, ((vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3)^2+
            ((vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3)^2 <= vmposmax^2)
    JuMP.@constraint(pm.model, ((vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3)^2+
            ((vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3)^2 >= vmposmin^2)

    # PMD.sol(pm, nw, :bus, bus_id)[:vmpossqr] = vrepos^2+vimpos^2
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_zero_sequence(pm::PMD.AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vmzeromax::Real)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmzerosqr)
        PMD.var(pm, PMD.nw_id_default)[:vmzerosqr] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMD.var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [PMD.var(pm, nw, :va, bus_id)[i] for i in 1:3]
    ### real and imaginary components of U+
    vrezero = JuMP.NonlinearExpr[]
    vimzero = JuMP.NonlinearExpr[]
    vmzerosqr = JuMP.NonlinearExpr[]
    vrezero = JuMP.@expression(pm.model,
        (vm_a*cos(va_a) + vm_b*cos(va_b) + vm_c*cos(va_c))/3
    )
    vimzero = JuMP.@expression(pm.model,
        (vm_a*sin(va_a) + vm_b*sin(va_b) + vm_c*sin(va_c))/3
    )

    ### square of magnitude of U+, |U+|^2
    vmzerosqr = JuMP.@expression(pm.model, vrezero^2+vimzero^2)
    ### finally, apply constraint
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
function constraint_mc_bus_voltage_magnitude_vuf(pm::PMD.AbstractUnbalancedACPModel, nw::Int, bus_id::Int, vufmax::Real)
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmpossqr)
        PMD.var(pm, PMD.nw_id_default)[:vmpossqr] = Dict{Int, Any}()
    end
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vmnegsqr)
        PMD.var(pm, PMD.nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
    end
    if !haskey(PMD.var(pm, PMD.nw_id_default), :vuf)
        PMD.var(pm, PMD.nw_id_default)[:vuf] = Dict{Int, Any}()
    end
    (vm_a, vm_b, vm_c) = [PMD.var(pm, nw, :vm, bus_id)[i] for i in 1:3]
    (va_a, va_b, va_c) = [PMD.var(pm, nw, :va, bus_id)[i] for i in 1:3]
    a = exp(im*2*pi/3)
    # real and imag functions cannot be used in NLexpressions, so precalculate
    are = real(a)
    aim = imag(a)
    a2re = real(a^2)
    a2im = imag(a^2)

    # real and imaginary components of U+
    vrepos = JuMP.NonlinearExpr[]
    vimpos = JuMP.NonlinearExpr[]
    vmpossqr = JuMP.NonlinearExpr[]
    vrepos = JuMP.@expression(pm.model,
        (vm_a*cos(va_a) + are*vm_b*cos(va_b) - aim*vm_b*sin(va_b) + a2re*vm_c*cos(va_c) - a2im*vm_c*sin(va_c))/3
    )
    vimpos = JuMP.@expression(pm.model,
        (vm_a*sin(va_a) + are*vm_b*sin(va_b) + aim*vm_b*cos(va_b) + a2re*vm_c*sin(va_c) + a2im*vm_c*cos(va_c))/3
    )
    # square of magnitude of U+, |U+|^2
    vmpossqr = JuMP.@expression(pm.model, vrepos^2+vimpos^2)

    # real and imaginary components of U-
    vreneg = JuMP.NonlinearExpr[]
    vimneg = JuMP.NonlinearExpr[]
    vmnegsqr = JuMP.NonlinearExpr[]
    vreneg = JuMP.@expression(pm.model,
        (vm_a*cos(va_a) + a2re*vm_b*cos(va_b) - a2im*vm_b*sin(va_b) + are*vm_c*cos(va_c) - aim*vm_c*sin(va_c))/3
    )
    vimneg = JuMP.@expression(pm.model,
        (vm_a*sin(va_a) + a2re*vm_b*sin(va_b) + a2im*vm_b*cos(va_b) + are*vm_c*sin(va_c) + aim*vm_c*cos(va_c))/3
    )
    # square of magnitude of U-, |U-|^2
    vmnegsqr = JuMP.@expression(pm.model, vreneg^2+vimneg^2)

    vuf = JuMP.NonlinearExpr[]
    vuf = JuMP.@expression(pm.model, vreneg / vimneg)
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmnegsqr <= vufmax^2*vmpossqr)
    # DEBUGGING: save references for post check
    #PMD.var(pm, PMD.nw_id_default, :vmpossqr)[bus_id] = vmpossqr
    #PMD.var(pm, PMD.nw_id_default, :vmnegsqr)[bus_id] = vmnegsqr

    PMD.sol(pm, nw, :bus, bus_id)[:vuf] = vuf
end


# """
# Bounds the current magnitude at both from and to side of a branch
# `cr[f_idx]^2 + ci[f_idx]^2 <= c_rating_a^2`
# `cr[t_idx]^2 + ci[t_idx]^2 <= c_rating_a^2`
# """
# function constraint_mc_current_limit(pm::PMD.AbstractUnbalancedACPModel, nw::Int, f_idx::Tuple{Int,Int,Int}, f_connections::Vector{Int}, t_connections::Vector{Int}, c_rating_a::Vector{<:Real})
#     (l, f_bus, t_bus) = f_idx
#     t_idx = (l, t_bus, f_bus)

#     crf =  [PMD.var(pm, nw, :cr, f_idx)[c] for c in f_connections]
#     cif =  [PMD.var(pm, nw, :ci, f_idx)[c] for c in f_connections]

#     crt =  [PMD.var(pm, nw, :cr, t_idx)[c] for c in t_connections]
#     cit =  [PMD.var(pm, nw, :ci, t_idx)[c] for c in t_connections]

#     JuMP.@constraint(pm.model, crf.^2 + cif.^2 .<= c_rating_a.^2)
#     JuMP.@constraint(pm.model, crt.^2 + cit.^2 .<= c_rating_a.^2)
# end
