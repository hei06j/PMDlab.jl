using Pkg
Pkg.activate("./")
using PMDlab
using CSV
using DataFrames


# define folder locations
pkg_dir = dirname(pathof(PMDlab))
cd( joinpath(pkg_dir, "..", "data") )

ENWL_ORIGINAL_DIR = "enwl_original"
ENWL_4W_EMBD_DIR = "enwl_4w"

LINE_MODELS_DIR = "urquhart_line_models"
LOAD_PROFILE_DIR = "Load_profiles"


mkpath("$LOAD_PROFILE_DIR")
embedded_load_index = 10 # To select the load timestep out of 288 (5-minutely)
embedded_load_profile = []
mult_factor = 6 # To increas the load values if necessary
load_profile = CSV.read("./load_profiles_original/Summer_Load_Profiles.csv", DataFrame, header=false)
n_profiles = length(names(load_profile))
for i in 1:4  # there is only 100 loads available, whereas some feeders require up to 400 load profiles, so we duplicate the existing ones
    for (h, col) in enumerate(names(load_profile))
        # output_file = "$LOAD_PROFILE_DIR/Load$(h)_5min.csv"
        output_file = "$LOAD_PROFILE_DIR/Load$(h+(i-1)*n_profiles)_5min.csv"
        df = DataFrame(col=load_profile[!,"$col"]) .* mult_factor
        open(output_file, "w") do f
            CSV.write(f, df, writeheader=false)
        end
        append!(embedded_load_profile, df[embedded_load_index,1])
    end
end

mutable struct embedded_load
    is_embedded::Bool
    profile::Vector{Float64}
end


# load linecode data and calculate positive sequence impedance
target_lcs_4w = Dict{Int,NamedTuple}()

for (k,fn) in enumerate(readdir("$LINE_MODELS_DIR/input_files"))
    if fn != "current_limits.md"
        RX = Matrix(CSV.read("$LINE_MODELS_DIR/input_files/$fn", DataFrame, header=false))
        Z = [RX[i, (j-1)*2+1]+im*RX[i, (j-1)*2+2] for i in 1:5, j in 1:5]
        Zfw = PMDlab.kron_reduce(Z, 5)
        Zkr = PMDlab.kron_reduce(Zfw, 4)
        Zseq = PMDlab.seq_transform_zpn(Zkr)
        zpos = Zseq[2,2]
        
        lc_dss_def = Vector{String}()
        push!(lc_dss_def, "! original name: $(fn[1:end-4])")
        push!(lc_dss_def, "New LineCode.LC$k nphases=$(size(Zfw)[1]) Units=km")
        push!(lc_dss_def, "~ Rmatrix="*PMDlab._mat_to_dss_string(real.(Zfw)))
        push!(lc_dss_def, "~ Xmatrix="*PMDlab._mat_to_dss_string(imag.(Zfw)))
        push!(lc_dss_def, "~ Cmatrix="*PMDlab._mat_to_dss_string(zeros(size(Zfw)...)))
        target_lcs_4w[k] = (zpos=zpos, zfw=Zfw, dss=lc_dss_def)

    end
end


# save file with all linecodes
dss_all_target_lcs_4w = PMDlab._combine_dss_subset_target_lcs(target_lcs_4w, 1:length(target_lcs_4w))
open("$LINE_MODELS_DIR/linecodes_4w.dss", "w") do f
    write(f, join(dss_all_target_lcs_4w, "\n"))
end

with_transformer = false

# iterate over all cases and generate 4w, phase_to_neutral, kron_reduced, and modified_phase_to_neutral cases
for network_dir in readdir(ENWL_ORIGINAL_DIR)
    mn = match(r"network_(\d+)$"i, network_dir)
    if !isnothing(mn)
        n = mn.captures[1]
        for feeder_dir in readdir("$ENWL_ORIGINAL_DIR/$network_dir")
            mf = match(r"feeder_(\d+)$"i, feeder_dir)
            if !isnothing(mf)
                f = mf.captures[1]
                case_dir = "$ENWL_ORIGINAL_DIR/$network_dir/$feeder_dir"

                # 4 wire, snap (embedded) load
                n_wire = "4"
                output_folder = "$ENWL_4W_EMBD_DIR//$network_dir//$feeder_dir"
                cct_name = network_dir*"_"*feeder_dir*"_4wire"
                mkpath(output_folder)
                PMDlab.build_case(case_dir, target_lcs_4w, output_folder, n_wire, cct_name, embedded_load(true,embedded_load_profile), with_transformer)
            end
        end
    end
end

cd("..")