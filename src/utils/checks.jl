
function check_active_bounds(result3w_ivr, math3w)
    vm_bounds = Dict()
    vmneg_bounds = Dict()
    vmzero_bounds = Dict()
    vmpos_bounds = Dict()
    vuf_bounds = Dict()
    for (i, bus) in result3w_ivr["solution"]["bus"]
        vm_bounds[i] = math3w["bus"][i]["vmin"] .< sqrt.(bus["vr"].^2 .+ bus["vi"].^2) .< math3w["bus"][i]["vmax"]
        vmneg_bounds[i] = haskey(bus, "vmnegsqr") ? sqrt(bus["vmnegsqr"]) < math3w["bus"][i]["vm_seq_neg_max"] : true
        vmzero_bounds[i] = haskey(bus, "vmzerosqr") ? sqrt(bus["vmzerosqr"]) < math3w["bus"][i]["vm_seq_zero_max"] : true
        vmpos_bounds[i] = haskey(bus, "vmpossqr") ? sqrt(bus["vmpossqr"]) < math3w["bus"][i]["vm_seq_pos_max"] : true
        vuf_bounds[i] = haskey(bus, "vmzerosqr") ? sqrt(bus["vmnegsqr"]/bus["vmpossqr"]) < math3w["bus"][i]["vm_vuf_max"] : true
        if !all(vm_bounds[i])
            print("bus $i with active $(vm_bounds[i]) \n vm=$(sqrt.(bus["vr"].^2 .+ bus["vi"].^2)) and vmin=$(math3w["bus"][i]["vmin"]), vmax=$(math3w["bus"][i]["vmax"]) \n\n")
        end
        if !all(vmneg_bounds[i])
            print("bus $i with active $(vmneg_bounds[i]) \n vmneg=$(sqrt(bus["vmnegsqr"])) and vmneg_max=$(math3w["bus"][i]["vm_seq_neg_max"]) \n\n")
        end
        if !all(vmzero_bounds[i])
            print("bus $i with active $(vmzero_bounds[i]) \n vmzero=$(sqrt(bus["vmzerosqr"])) and vmzero_max=$(math3w["bus"][i]["vm_seq_zero_max"]) \n\n")
        end
        if !all(vmpos_bounds[i])
            print("bus $i with active $(vmpos_bounds[i]) \n vmpos=$(sqrt(bus["vmpossqr"])) and vmpos_max=$(math3w["bus"][i]["vm_seq_pos_max"]) \n\n")
        end
        if !all(vuf_bounds[i])
            print("bus $i with active $(vuf_bounds[i]) \n vuf=$(sqrt(bus["vmnegsqr"]/bus["vmpossqr"])) and vuf_max=$(math3w["bus"][i]["vm_vuf_max"]) \n\n")
        end
    end


    cm_bounds = Dict()
    for (i, branch) in result3w_ivr["solution"]["branch"]
        cm_bounds[i] = sqrt.(branch["cr_fr"].^2 .+ branch["ci_fr"].^2) .<= math3w["branch"][i]["c_rating_a"]
        if !all(cm_bounds[i])
            print("branch $i with active $(cm_bounds[i]) \n c_from=$(sqrt.(branch["cr_fr"].^2 .+ branch["ci_fr"].^2)) and c_rating_a=$(math3w["branch"][i]["c_rating_a"]) \n\n")
        end
    end


    pg_bounds = Dict()
    qg_bounds = Dict()
    for (i, gen) in result3w_ivr["solution"]["gen"]
        pg_bounds[i] = math3w["gen"][i]["pmin"] .< gen["pg"] .< math3w["gen"][i]["pmax"]
        qg_bounds[i] = math3w["gen"][i]["qmin"] .< gen["qg"] .< math3w["gen"][i]["qmax"]
        if !all(pg_bounds[i])
            print("gen $i with active $(pg_bounds[i]) \n pg=$(gen["pg"]) with pmin=$(math3w["gen"][i]["pmin"]) and pmax=$(math3w["gen"][i]["pmax"]) \n\n")
        end
        if !all(qg_bounds[i])
            print("gen $i with active $(qg_bounds[i]) \n qg=$(gen["qg"]) with qmin=$(math3w["gen"][i]["qmin"]) and qmax=$(math3w["gen"][i]["qmax"]) \n\n")
        end
    end
end
