using Pkg
Pkg.activate("./")
using PMDlab
using PowerModelsDistribution
using Ipopt
using JuMP
using LinearAlgebra

# optimizer = Ipopt.Optimizer
optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0, "sb"=>"yes","warm_start_init_point"=>"yes", "max_iter"=>100)
optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "max_iter"=>100)

const PMD = PowerModelsDistribution

## For all the test cases, the following new constraints are checked:
### - negative sequence bound should be binding in some of the testcase
### - current limits also should be binding in some of the testcase
ENWL_3W_EMBD_DIR = "data/three-wire"

vneg_vio_acp = []
cm_vio_acp = []

vneg_vio_acr = []
cm_vio_acr = []

vneg_vio_ivr = []
cm_vio_ivr = []

using Plots
plt = plot()
for network_dir in readdir(ENWL_3W_EMBD_DIR)
    mn = match(r"network_(\d+)$"i, network_dir)
    if !isnothing(mn)
        n = mn.captures[1]
        for feeder_dir in readdir("$ENWL_3W_EMBD_DIR/$network_dir")
            mf = match(r"feeder_(\d+)$"i, feeder_dir)
            if !isnothing(mf)
                f = mf.captures[1]
                case_dir = "$ENWL_3W_EMBD_DIR/$network_dir/$feeder_dir"
                @show case_dir

                eng3w = parse_file(case_dir*"/Master.dss", transformations=[transform_loops!])
                PMDlab.augment_eng_3wire!(eng3w; line_current_rating=true, sbase=1)

                math3w = transform_data_model(eng3w, kron_reduce=true, phase_project=true)
                PMDlab.augment_math_3wire!(math3w; relax_vsource_vm=true, Vsequence_bounds=true, cost_multiplier=1000)  # chaning some of the input data

                result3w_acp = PMD.solve_mc_model(math3w, ACPUPowerModel, optimizer, PMDlab.build_mc_opf)
                vneg_acp = [math3w["bus"][i]["vm_seq_neg_max"]-sqrt(JuMP.value(bus["vmnegsqr"])+1E-10) for (i,bus) in result3w_acp["solution"]["bus"] if haskey(bus, "vmnegsqr")]
                cm_acp = [math3w["branch"][i]["c_rating_a"] .- sqrt.((branch["pf"].^2+branch["qf"].^2)./result3w_acp["solution"]["bus"]["$(math3w["branch"][i]["f_bus"])"]["vm"].^2) for (i,branch) in result3w_acp["solution"]["branch"]]
                append!(vneg_vio_acp, minimum(vneg_acp))
                append!(cm_vio_acp, minimum(cm_acp))

                add_start_voltage!(math3w, coordinates=:rectangular, explicit_neutral=false)
                result3w_acr = PMD.solve_mc_model(math3w, ACRUPowerModel, optimizer, PMDlab.build_mc_opf)
                vneg_acr = [math3w["bus"][i]["vm_seq_neg_max"]-sqrt(bus["vmnegsqr"]+1E-10) for (i,bus) in result3w_acr["solution"]["bus"] if haskey(bus, "vmnegsqr")]
                cm_acr = [math3w["branch"][i]["c_rating_a"] .- sqrt.((branch["pf"].^2+branch["qf"].^2)./(result3w_acr["solution"]["bus"]["$(math3w["branch"][i]["f_bus"])"]["vr"].^2+result3w_acr["solution"]["bus"]["$(math3w["branch"][i]["f_bus"])"]["vi"].^2)) for (i,branch) in result3w_acr["solution"]["branch"]]
                append!(vneg_vio_acr, minimum(vneg_acr))
                append!(cm_vio_acr, minimum(cm_acr))

                add_start_voltage!(math3w, coordinates=:rectangular, explicit_neutral=false)
                result3w_ivr = PMD.solve_mc_model(math3w, IVRUPowerModel, optimizer, PMDlab.build_mc_opf)
                vneg_ivr = [math3w["bus"][i]["vm_seq_neg_max"]-sqrt(bus["vmnegsqr"]+1E-10) for (i,bus) in result3w_ivr["solution"]["bus"] if haskey(bus, "vmnegsqr")]
                cm_ivr = [math3w["branch"][i]["c_rating_a"] .- abs.(branch["cr_fr"]+im*branch["ci_fr"]) for (i,branch) in result3w_ivr["solution"]["branch"]]
                append!(vneg_vio_ivr, minimum(vneg_ivr))
                append!(cm_vio_ivr, minimum(cm_ivr))

            end
        end
    end
end

cm_vio = scatter(cm_vio_ivr, label="IVR", title="Minimum (Crating - Cm) for each network")
scatter!(cm_vio_acr, label="ACR")
scatter!(cm_vio_acp, label="ACP")
savefig(cm_vio, "Figures/Current_limits_violation_validation.pdf")

vneg_vio = scatter(vneg_vio_ivr, label="IVR", title="Minimum (Vneg_max - Vneg) for each network")
scatter!(vneg_vio_acr, label="ACR")
scatter!(vneg_vio_acp, label="ACP")
savefig(vneg_vio, "Figures/Vneg_limits_violation_validation.pdf")


## Inspect sequence voltage components

vpos_acp = [sqrt(JuMP.value(bus["vmpossqr"])) for (i,bus) in result3w_acp["solution"]["bus"] if haskey(bus, "vmpossqr")]
vneg_acp = [sqrt(JuMP.value(bus["vmnegsqr"])) for (i,bus) in result3w_acp["solution"]["bus"] if haskey(bus, "vmnegsqr")]
vzero_acp = [sqrt(JuMP.value(bus["vmzerosqr"]+1E-10)) for (i,bus) in result3w_acp["solution"]["bus"] if haskey(bus, "vmzerosqr")]
vuf_acp = [abs(JuMP.value(bus["vuf"]))*100 for (i,bus) in result3w_acp["solution"]["bus"] if haskey(bus, "vuf")]

vpos_acr = [sqrt(bus["vmpossqr"]) for (i,bus) in result3w_acr["solution"]["bus"] if haskey(bus, "vmpossqr")]
vneg_acr = [sqrt(bus["vmnegsqr"]) for (i,bus) in result3w_acr["solution"]["bus"] if haskey(bus, "vmnegsqr")]
vzero_acr = [sqrt(bus["vmzerosqr"]+1E-10) for (i,bus) in result3w_acr["solution"]["bus"] if haskey(bus, "vmzerosqr")]
vuf_acr = [abs(JuMP.value(bus["vuf"]))*100 for (i,bus) in result3w_acr["solution"]["bus"] if haskey(bus, "vuf")]

vpos_ivr = [sqrt(bus["vmpossqr"]) for (i,bus) in result3w_ivr["solution"]["bus"] if haskey(bus, "vmpossqr")]
vneg_ivr = [sqrt(bus["vmnegsqr"]) for (i,bus) in result3w_ivr["solution"]["bus"] if haskey(bus, "vmnegsqr")]
vzero_ivr = [sqrt(bus["vmzerosqr"]+1E-10) for (i,bus) in result3w_ivr["solution"]["bus"] if haskey(bus, "vmzerosqr")]
vuf_ivr = [abs(JuMP.value(bus["vuf"]))*100 for (i,bus) in result3w_ivr["solution"]["bus"] if haskey(bus, "vuf")]

using Plots

scatter(vpos_acp)
scatter!(vpos_acr)
scatter!(vpos_ivr)

scatter(vneg_acp)
scatter!(vneg_acr)
scatter!(vneg_ivr)

scatter(vzero_acp)
scatter!(vzero_acr)
scatter!(vzero_ivr)

# scatter(vuf_acp)
# scatter!(vuf_acr)
# scatter!(vuf_ivr)


## Inspect current limits

cm_ivr = [math3w["branch"][i]["c_rating_a"] .- abs.(branch["cr_fr"]+im*branch["ci_fr"]) for (i,branch) in result3w_ivr["solution"]["branch"]]
cm_acr = [math3w["branch"][i]["c_rating_a"] .- sqrt.((branch["pf"].^2+branch["qf"].^2)./(result3w_acr["solution"]["bus"]["$(math3w["branch"][i]["f_bus"])"]["vr"].^2+result3w_acr["solution"]["bus"]["$(math3w["branch"][i]["f_bus"])"]["vi"].^2)) for (i,branch) in result3w_acr["solution"]["branch"]]
cm_acp = [math3w["branch"][i]["c_rating_a"] .- sqrt.((branch["pf"].^2+branch["qf"].^2)./result3w_acp["solution"]["bus"]["$(math3w["branch"][i]["f_bus"])"]["vm"].^2) for (i,branch) in result3w_acp["solution"]["branch"]]

cm_ivr = [branch["c_rating_a"] for (i,branch) in math3w["branch"]]

minimum(cm_ivr)
minimum(cm_acr)
minimum(cm_acp)

scatter(cm_ivr, label=false)
scatter(cm_acr, label=false)
scatter(cm_acp, label=false)

## Inspect voltage magnitudes
# [result3w_acp["solution"]["gen"]["1"]["pg"]  result3w_acr["solution"]["gen"]["1"]["pg"] result3w_ivr["solution"]["gen"]["1"]["pg"]]
# [result3w_acp["solution"]["gen"]["2"]["pg"]  result3w_acr["solution"]["gen"]["2"]["pg"] result3w_ivr["solution"]["gen"]["2"]["pg"]]

# result3w_acp["solution"]["bus"]["$(ref_bus[1])"]
# abs.(result3w_acr["solution"]["bus"]["$(ref_bus[1])"]["vr"] .+ im* result3w_acr["solution"]["bus"]["$(ref_bus[1])"]["vi"])
# angle.(result3w_acr["solution"]["bus"]["$(ref_bus[1])"]["vr"] .+ im* result3w_acr["solution"]["bus"]["$(ref_bus[1])"]["vi"])

# abs.(result3w_ivr["solution"]["bus"]["$(ref_bus[1])"]["vr"] .+ im* result3w_ivr["solution"]["bus"]["$(ref_bus[1])"]["vi"])
# angle.(result3w_ivr["solution"]["bus"]["$(ref_bus[1])"]["vr"] .+ im* result3w_ivr["solution"]["bus"]["$(ref_bus[1])"]["vi"])

# using Plots

# vma_acr = [abs.(result3w_acr["solution"]["bus"][i]["vr"][1] .+ im* result3w_acr["solution"]["bus"][i]["vi"][1]) for (i,bus) in math3w["bus"]]
# vmb_acr = [abs.(result3w_acr["solution"]["bus"][i]["vr"][2] .+ im* result3w_acr["solution"]["bus"][i]["vi"][2]) for (i,bus) in math3w["bus"]]
# vmc_acr = [abs.(result3w_acr["solution"]["bus"][i]["vr"][3] .+ im* result3w_acr["solution"]["bus"][i]["vi"][3]) for (i,bus) in math3w["bus"]]

# vma_acp = [result3w_acp["solution"]["bus"][i]["vm"][1] for (i,bus) in math3w["bus"]]
# vmb_acp = [result3w_acp["solution"]["bus"][i]["vm"][2] for (i,bus) in math3w["bus"]]
# vmc_acp = [result3w_acp["solution"]["bus"][i]["vm"][3] for (i,bus) in math3w["bus"]]

# vma_ivr = [abs.(result3w_ivr["solution"]["bus"][i]["vr"][1] .+ im* result3w_ivr["solution"]["bus"][i]["vi"][1]) for (i,bus) in math3w["bus"]]
# vmb_ivr = [abs.(result3w_ivr["solution"]["bus"][i]["vr"][2] .+ im* result3w_ivr["solution"]["bus"][i]["vi"][2]) for (i,bus) in math3w["bus"]]
# vmc_ivr = [abs.(result3w_ivr["solution"]["bus"][i]["vr"][3] .+ im* result3w_ivr["solution"]["bus"][i]["vi"][3]) for (i,bus) in math3w["bus"]]


# scatter(vma_acr)
# scatter!(vma_acp)
# scatter!(vma_ivr)

# scatter(vmb_acr)
# scatter!(vmb_acp)
# scatter!(vmb_ivr)

# scatter(vmc_acr)
# scatter!(vmc_acp)
# scatter!(vmc_ivr)
