"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_negative_sequence(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vmnegmax::Real)
    if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
        PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
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
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_positive_sequence(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vmposmax::Real)
    if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
        PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
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
end


"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_zero_sequence(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vmzeromax::Real)
    if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
        PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
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
end



"""
a = exp(im*2π/3)
U+ = (1*Ua + a*Ub a^2*Uc)/3
U- = (1*Ua + a^2*Ub a*Uc)/3
vuf = |U-|/|U+|
|U-| <= vufmax*|U+|
|U-|^2 <= vufmax^2*|U+|^2
"""
function constraint_mc_bus_voltage_magnitude_vuf(pm::PMD.AbstractUnbalancedIVRModel, nw::Int, bus_id::Int, vufmax::Real)
    if !haskey(PMD.var(pm, nw_id_default), :vmpossqr)
        PMD.var(pm, nw_id_default)[:vmpossqr] = Dict{Int, Any}()
        PMD.var(pm, nw_id_default)[:vmnegsqr] = Dict{Int, Any}()
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
    # finally, apply constraint
    JuMP.@constraint(pm.model, vmnegsqr <= vufmax^2*vmpossqr)
    # DEBUGGING: save references for post check
    #PMD.var(pm, nw_id_default, :vmpossqr)[bus_id] = vmpossqr
    #PMD.var(pm, nw_id_default, :vmnegsqr)[bus_id] = vmnegsqr
end
