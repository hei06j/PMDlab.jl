
function check_active_bounds(result3w, math3w)
    df = DataFrame(network=math3w["name"], solution=result3w["termination_status"],
        vm=false, va_pp=false, va_ij=false, 
        vmneg=false, vmzero=false, vmpos=false, vuf=false, 
        cm=false, pg=false, qg=false)
    
    vm_bounds = Dict()
    vmneg_bounds = Dict()
    vmzero_bounds = Dict()
    vmpos_bounds = Dict()
    vuf_bounds = Dict()
    va_pp_bounds = Dict()
    for (i, bus) in result3w["solution"]["bus"]
        vm_bounds[i] = math3w["bus"][i]["vmin"] .< sqrt.(bus["vr"].^2 .+ bus["vi"].^2) .< math3w["bus"][i]["vmax"]
        vmneg_bounds[i] = haskey(bus, "vmnegsqr") ? sqrt(bus["vmnegsqr"]) < math3w["bus"][i]["vm_seq_neg_max"] : true
        vmzero_bounds[i] = haskey(bus, "vmzerosqr") ? sqrt(bus["vmzerosqr"]) < math3w["bus"][i]["vm_seq_zero_max"] : true
        vmpos_bounds[i] = haskey(bus, "vmpossqr") ? sqrt(bus["vmpossqr"]) < math3w["bus"][i]["vm_seq_pos_max"] : true
        vuf_bounds[i] = haskey(bus, "vmzerosqr") ? sqrt(bus["vmnegsqr"]/bus["vmpossqr"]) < math3w["bus"][i]["vm_vuf_max"] : true

        T = [1 -1 0 ; 0 1 -1 ; -1 0 1]
        va_pp_bounds[i] = haskey(math3w["bus"][i], "va_delta") ? -math3w["bus"][i]["va_delta"] .< PMD._wrap_to_180(T * angle.(bus["vr"] .+ im * bus["vi"]) .* 180/pi .- 120)  .< math3w["bus"][i]["va_delta"] : true

        if !all(vm_bounds[i])
            df[1,"vm"] = true
            print("bus $i with active $(vm_bounds[i]) \n vm=$(sqrt.(bus["vr"].^2 .+ bus["vi"].^2)) and vmin=$(math3w["bus"][i]["vmin"]), vmax=$(math3w["bus"][i]["vmax"]) \n\n")
        end
        if !all(vmneg_bounds[i])
            df[1,"vmneg"] = true
            print("bus $i with active $(vmneg_bounds[i]) \n vmneg=$(sqrt(bus["vmnegsqr"])) and vmneg_max=$(math3w["bus"][i]["vm_seq_neg_max"]) \n\n")
        end
        if !all(vmzero_bounds[i])
            df[1,"vmzero"] = true
            print("bus $i with active $(vmzero_bounds[i]) \n vmzero=$(sqrt(bus["vmzerosqr"])) and vmzero_max=$(math3w["bus"][i]["vm_seq_zero_max"]) \n\n")
        end
        if !all(vmpos_bounds[i])
            df[1,"vmpos"] = true
            print("bus $i with active $(vmpos_bounds[i]) \n vmpos=$(sqrt(bus["vmpossqr"])) and vmpos_max=$(math3w["bus"][i]["vm_seq_pos_max"]) \n\n")
        end
        if !all(vuf_bounds[i])
            df[1,"vuf"] = true
            print("bus $i with active $(vuf_bounds[i]) \n vuf=$(sqrt(bus["vmnegsqr"]/bus["vmpossqr"])) and vuf_max=$(math3w["bus"][i]["vm_vuf_max"]) \n\n")
        end
        if !all(va_pp_bounds[i])
            df[1,"va_pp"] = true
            print("bus $i with active $(va_pp_bounds[i]) \n va_pp=$(PMD._wrap_to_180(T * angle.(bus["vr"] .+ im * bus["vi"]) .* 180/pi .- 120)) and va_delta=$(math3w["bus"][i]["va_delta"]) \n\n")
        end
    end
    
    va_ij_bounds = Dict()
    cm_bounds = Dict()
    for (i, branch) in math3w["branch"]
        cm_bounds[i] = sqrt.(result3w["solution"]["branch"][i]["cr_fr"].^2 .+ result3w["solution"]["branch"][i]["ci_fr"].^2) .<= branch["c_rating_a"]
        fbus = branch["f_bus"]
        tbus = branch["t_bus"]
        fbus_va = angle.(result3w["solution"]["bus"]["$fbus"]["vr"] .+ im * result3w["solution"]["bus"]["$fbus"]["vi"])
        tbus_va = angle.(result3w["solution"]["bus"]["$tbus"]["vr"] .+ im * result3w["solution"]["bus"]["$tbus"]["vi"])
        va_ij_bounds[i] = haskey(branch, "angmin") ? branch["angmin"] .< (fbus_va .- tbus_va) .< branch["angmax"] : true

        if !all(cm_bounds[i])
            df[1,"cm"] = true
            print("branch $i with active $(cm_bounds[i]) \n c_from=$(sqrt.(result3w["solution"]["branch"][i]["cr_fr"].^2 .+ result3w["solution"]["branch"][i]["ci_fr"].^2)) and c_rating_a=$(branch["c_rating_a"]) \n\n")
        end
        if !all(va_ij_bounds[i])
            df[1,"va_ij"] = true
            print("branch $i with active $(va_ij_bounds[i]) \n va_ij=$((fbus_va .- tbus_va) .* 180/pi) and angmin=$(branch["angmin"].*180/pi) and angmax=$(branch["angmax"].*180/pi) \n\n")
        end
    end


    pg_bounds = Dict()
    qg_bounds = Dict()
    for (i, gen) in result3w["solution"]["gen"]
        pg_bounds[i] = math3w["gen"][i]["pmin"] .< gen["pg"] .< math3w["gen"][i]["pmax"]
        qg_bounds[i] = math3w["gen"][i]["qmin"] .< gen["qg"] .< math3w["gen"][i]["qmax"]
        if !all(pg_bounds[i])
            df[1,"pg"] = true
            print("gen $i with active $(pg_bounds[i]) \n pg=$(gen["pg"]) with pmin=$(math3w["gen"][i]["pmin"]) and pmax=$(math3w["gen"][i]["pmax"]) \n\n")
        end
        if !all(qg_bounds[i])
            df[1,"qg"] = true
            print("gen $i with active $(qg_bounds[i]) \n qg=$(gen["qg"]) with qmin=$(math3w["gen"][i]["qmin"]) and qmax=$(math3w["gen"][i]["qmax"]) \n\n")
        end
    end

    check_dict = Dict("vm" => vm_bounds, "va_pp_bounds" => vm_bounds, "va_ij" => va_ij_bounds, 
                      "vmneg" => vmneg_bounds, "vmzero" => vmzero_bounds, "vmpos" => vmpos_bounds, "vuf" => vuf_bounds,
                      "cm" => cm_bounds, "pg" => pg_bounds, "qg" => qg_bounds)

    return df, check_dict
end

## have upper and lower bounds for checks