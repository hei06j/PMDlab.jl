_mat_to_dss_string(mat::Matrix{<:Real}) = "("*join([join(mat[i,:], " ") for i in 1:size(mat)[1]], " | ")*")"


"transforms a matrix (Z/Y) to phase-neutral impedance model"
function three_wire_reduce(Z::Matrix)
    @assert(all(size(Z).==4), "Matrix has to be 4x4!")
    A = [1 0 0 -1 ; 0 1 0 -1 ; 0 0 1 -1]
    Zpn = A * Z * A'
    return Zpn
end


"transforms a matrix (Z/Y) to modified phase-neutral impedance model"
function three_wire_modified_reduce(Z::Matrix)
    @assert(all(size(Z).==4), "Matrix has to be 4x4!")
    A = [1 0 0 -1 ; 0 1 0 -1 ; 0 0 1 -1]
    ZMpn = A * diagm(diag(Z)) * A'
    return ZMpn
end


"transforms a matrix (Z/Y) to kron-reduced impedance model"
function kron_reduce(Z::Matrix, n::Int)
    N = size(Z)[1]; @assert size(Z)[2]==N
    P = setdiff(1:N,  n)
    Zkr = Z[P,P]-(1/Z[n,n])*Z[P,[n]]*Z[[n],P]
    return Zkr
end


"transforms a matrix (Z/Y) to a sequence component basis"
function seq_transform_zpn(M)
    @assert(all(size(M).==3), "Matrix has to be 3x3!")
    α = exp(im*2*pi/3)
    A = (1/sqrt(3))*[1 1 1; 1 α^2 α; 1 α α^2]
    Ainv = (1/sqrt(3))*[1 1 1; 1 α α^2; 1 α^2 α]
    Mseq = Ainv*M*A
    return Mseq
end


function _find_target_lc(target_lcs, zpos)
    min_dist = Inf
    target = missing
    for (lc, lc_data) in target_lcs
        dist = abs(zpos-lc_data.zpos)
        if dist<min_dist
            target = lc
            min_dist = dist
        end
    end
    return target
end


function _add_lines_with_header!(lines, add_lines, header)
    # insert a blank line if not present
    if !isempty(lines)
        if !isempty(lines[end])
            push!(lines, "")
        end
    end

    width = maximum(length.(add_lines))
    push!(lines, "! "*header)
    push!(lines, "! "*repeat("-", width-2))
    append!(lines, add_lines)
    push!(lines, "! "*repeat("-", width-2))
end


function _combine_dss_subset_target_lcs(target_lcs::Dict{Int,NamedTuple}, subset)
    lines = Vector{String}()
    for s in subset
        append!(lines, target_lcs[s].dss)
        push!(lines, "")
    end
    return lines[1:end-1] # remove final newline
end


function build_master_dss(output_folder, cct_name, n_wire, embedded_load, with_transformer)
    master_file = "$output_folder/Master.dss"

    # build the case as a vector of strings, each representing a line
    lines_master = Vector{String}()
    
    # license information
    add_license_text!(lines_master)

    # add the circuit header
    push!(lines_master, "Clear")
    push!(lines_master, "Set DefaultBaseFreq=50")

    # add the circuit elements

    if with_transformer
    # push!(lines_master, "New Circuit.ENWL_$cct_name BasekV=11 pu=1.00 ISC3=3000 ISC1=2500")
        push!(lines_master, "New Circuit.ENWL_$cct_name BasekV=11 pu=1.00 ISC3=100000 ISC1=100000")
        push!(lines_master, "")
        push!(lines_master, "Redirect Transformers.txt")
    else
        push!(lines_master, "New Circuit.ENWL_$cct_name BasekV=0.416 pu=1.00 ISC3=100000 ISC1=100000")
        if n_wire == "4"
            push!(lines_master, "New Reactor.Grounding phases=1 bus1=sourcebus.4 bus2=sourcebus.0 R=0.0 X=1E-10")
        end
        push!(lines_master, "")
    end

    push!(lines_master, "Redirect LineCode.txt")
    push!(lines_master, "Redirect Lines.txt")
    if !embedded_load.is_embedded                # Networks with snap load do not need loadshape.txt and monitors.txt
        push!(lines_master, "Redirect LoadShapes.txt")
        push!(lines_master, "Redirect Monitors.txt")
    end
    push!(lines_master, "Redirect Loads.txt")
    
    # add solve 
    push!(lines_master, "")
    push!(lines_master, "New Energymeter.substation Element=Line.LINE1 1")
    push!(lines_master, "")
    if embedded_load.is_embedded
        push!(lines_master, "Set mode=Snap")
    else
        push!(lines_master, "Set mode=Daily number=288 stepsize=5m")
    end
    push!(lines_master, "Solve")
    push!(lines_master, "Closedi")
    
    # write the lines to the output file
    open(master_file, "w") do f
        write(f, join(lines_master, "\n"))
    end
    return nothing
end


function add_license_text!(lines_master)
    # license information
    push!(lines_master, "! Original network data and time series from ENWL - Low Voltage Network Solutions project (in OpenDSS format)")
    push!(lines_master, "!   https://www.enwl.co.uk/future-energy/innovation/smaller-projects/low-carbon-networks-fund/low-voltage-network-solutions/")
    push!(lines_master, "")
    push!(lines_master, "! Adapted with length-normalized four-wire impedance data from: ")
    push!(lines_master, "!   Urquhart, Andrew J., and Murray Thomson. 2019. “Cable Impedance Data”. figshare. https://hdl.handle.net/2134/15544.")
    push!(lines_master, "!   Creative Commons Attribution-NonCommercial 3.0 Unported License." )
    push!(lines_master, "")
    push!(lines_master, "! Adaptation process described in  ")
    push!(lines_master, "!   'Distribution Network Modeling: From Simulation Towards Optimization, Sander Claeys, Phd Dissertation, KU Leuven, Belgium 2021")
    push!(lines_master, "")
    push!(lines_master, "! Impedance transformations described in " )
    push!(lines_master, "!   Frederik Geth, Rahmat Heidari, and Arpan Koirala. 2022. Computational analysis of impedance transformations for four-wire power networks with sparse neutral grounding. In Proceedings of the Thirteenth ACM International Conference on Future Energy Systems (e-Energy '22). Association for Computing Machinery, New York, NY, USA, 105–113. https://doi.org/10.1145/3538637.3538844

! This data is otherwise licensed under [CC BY 4.0  Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/deed.en).
! The data can be retrieved from https://data.csiro.au/
" )
    push!(lines_master, "")
    return nothing
end


# add the transformer definition
function build_transformers_txt(case_dir, output_folder, n_wire)
    if n_wire ∉ ["1", "3", "4"]
        throw("Number of wires is not valid")
    end

    transformer_file = "$output_folder/Transformers.txt"
    lines_tr = Vector{String}()
    if n_wire == "4"
        trans_defs = ["New Transformer.TR1 Buses=[SourceBus.1.2.3 1.1.2.3.4] Conns=[Delta Wye] kVs=[11 0.416] kVAs=[800 800] XHL=1 phases=3"]
        push!(trans_defs, "New Reactor.Grounding phases=1 bus1=1.4 bus2=1.0 R=0.0 X=1E-10")
    elseif n_wire == "3"
        trans_defs = ["New Transformer.TR1 Buses=[SourceBus.1.2.3 1.1.2.3.0] Conns=[Delta Wye] kVs=[11 0.416] kVAs=[800 800] XHL=1 phases=3"]
    elseif n_wire == "1"
        trans_defs = ["New Transformer.TR1 Buses=[SourceBus 1] Conns=[Delta Wye] kVs=[11 0.416] kVAs=[800 800] XHL=1 phases=1"]
    end
    _add_lines_with_header!(lines_tr, trans_defs, "Transformer definition")
    
    open(transformer_file, "w") do f
        write(f, join(lines_tr, "\n"))
    end
    return nothing
end


# load the old linecode definitions
function build_linecode_txt(case_dir, target_lcs, output_folder)
    linecode_file = "$output_folder/LineCode.txt"

    lines_lc = Vector{String}()
    lc_old_defs = open("$case_dir/LineCode.txt", "r") do f
        readlines(f)
    end

    # find the target new linecode for each old linecode
    lc_mapping = Dict()
    for lc_old_def in lc_old_defs
        lc, rpos_str, xpos_str = match(r"^New LineCode.([^ ]+).*R1=([^ ]+).*X1=([^ ]+)", lc_old_def).captures
        # get the positive sequence impedance
        zpos = parse(Float64, rpos_str) + im*parse(Float64, xpos_str)
        # find the target new linecode which is has the 'closest' positive sequence impedance
        lc_mapping[lc] = _find_target_lc(target_lcs, zpos)
    end

    # add the relevant target new linecode definitions
    lc_new_defs = _combine_dss_subset_target_lcs(target_lcs, sort(unique(values(lc_mapping))))
    _add_lines_with_header!(lines_lc, lc_new_defs, "Subset of available target linecodes, which are closest to the original ones in terms of positive sequence impedance")
    
    open(linecode_file, "w") do f
        write(f, join(lines_lc, "\n"))
    end

    return lc_mapping
end


# load the line definitions, and modify them to refer to the applicable target new linecode
function build_lines_txt(case_dir, output_folder, target_lcs, lc_mapping, n_wire, with_transformer)
    if n_wire ∉ ["1", "3", "4"]
        throw("Number of wires is not valid")
    end

    lines_file = "$output_folder/Lines.txt"
    lines_line = Vector{String}()
    line_defs = open("$case_dir/Lines.txt", "r") do f
        readlines(f)
    end

    if !with_transformer
        prefix, bus1, bus2, nph, lc, suffix = match(r"(.*)Bus1=(\d+) Bus2=(\d+) phases=(\d) Linecode=([^ ]+)(.*)", line_defs[1]).captures
        line_defs = ["New Line.LINE0 Bus1=sourcebus Bus2=1 phases=3 Linecode=$(lc) Length=1 Units=m"; line_defs]
    end

    for k in 1:length(line_defs)
        line_def = line_defs[k]
        prefix, bus1, bus2, nph, lc, suffix = match(r"(.*)Bus1=(.*) Bus2=(\d+) phases=(\d) Linecode=([^ ]+)(.*)", line_def).captures
        lc_target = lc_mapping[lc]
        if n_wire == "4"
            line_defs[k] = "$(prefix)Bus1=$bus1.1.2.3.4 Bus2=$bus2.1.2.3.4 phases=4 Linecode=LC$(lc_target)$(suffix)"
        elseif n_wire == "3"
            line_defs[k] = "$(prefix)Bus1=$bus1.1.2.3 Bus2=$bus2.1.2.3 phases=3 Linecode=LC$(lc_target)$(suffix)"
        elseif n_wire == "1"
            line_defs[k] = "$(prefix)Bus1=$bus1 Bus2=$bus2 phases=1 Linecode=LC$(lc_target)$(suffix)"
        end
    end
    _add_lines_with_header!(lines_line, line_defs, "Imported from Lines.txt, with modified linecodes")
    
    open(lines_file, "w") do f
        write(f, join(lines_line, "\n"))
    end
    return nothing
end


# load the loadshape definitions, and modify them to point to correct folder
function build_loadshapes_txt(case_dir, output_folder)
    loadshapes_file = "$output_folder/LoadShapes.txt"
    lines_loadshape = Vector{String}()
    loadshape_defs = open("$case_dir/LoadShapes.txt", "r") do f
        readlines(f)
    end
    for k in 1:length(loadshape_defs)
        loadshape_def = loadshape_defs[k]
        prefix, dir, suffix = match(r"(.*)C:(.*)profiles\\(.*)", loadshape_def).captures
        loadshape_defs[k] = "$(prefix)../../../Load_profiles/$(suffix)"
    end
    _add_lines_with_header!(lines_loadshape, loadshape_defs, "Imported from LoadShapes.txt, with modified directory")
    
    open(loadshapes_file, "w") do f
        write(f, join(lines_loadshape, "\n"))
    end

    return nothing
end


# load the load definitions, and modify them to have the fourth wire as a return terminal instead of ground
function build_loads_txt(case_dir, output_folder, n_wire, embedded_load)
    if n_wire ∉ ["1", "3", "4"]
        throw("Number of wires is not valid")
    end

    loads_file = "$output_folder/Loads.txt"
    lines_load = Vector{String}()
    load_defs = open("$case_dir/Loads.txt", "r") do f
        readlines(f)
    end
    for k in 1:length(load_defs)
        load_def = load_defs[k]
        prefix, nph, bus_def, kV, kW, suffix = match(r"(.*)Phases=(\d) Bus1=([^ ]+) kV=([^ ]+) kW=([^ ]+)(.*)", load_def).captures
        if embedded_load.is_embedded
            suffix, Daily =  match(r"(.*) Daily=Shape_(.*)", suffix).captures
            kW = embedded_load.profile[parse(Int64, Daily), 1]
            # kW = CSV.read("Load_profiles/Load$(Daily)_5min.csv", DataFrame, header=false)[load_index,"Column1"]
        end

        if n_wire == "4"
            if nph=="3"
                load_defs[k] = "$(prefix)Phases=$nph Bus1=$(bus_def).1.2.3.4 kV=$kV kW=$kW $(suffix)"
            elseif nph=="1"
                bus, phase = split(bus_def, ".")
                load_defs[k] = "$(prefix)Phases=$nph Bus1=$(bus).$(phase).4 kV=$kV kW=$kW $(suffix)"
            else
                error("Encountered load with $nph phases.")
            end
        elseif n_wire == "3"
            if nph == "3"
                load_defs[k] = "$(prefix)Phases=$nph Bus1=$(bus_def).1.2.3 kV=$kV kW=$kW $(suffix)"
            elseif nph == "1"
                bus, phase = split(bus_def, ".")
                load_defs[k] = "$(prefix)Phases=$nph Bus1=$(bus).$(phase) kV=$kV kW=$kW $(suffix)"
            else
                error("Encountered load with $nph phases.")
            end
        elseif n_wire == "1"
            if nph == "3"
                kV = parse(Float64, kV)/sqrt(3)
                # @show kW, typeof(kW), String(kW), typeof(String(kW))
                kW = typeof(kW) == Float64 ? kW/3 : parse(Float64, kW)/3
                # kW = parse(Float64, String(kW))/3 
                # kv, kw, suffix3 = match(r" kV=([^ ]+) kW=(\d)(.*)", suffix).captures
                load_defs[k] = "$(prefix)Phases=1 Bus1=$(bus_def) kV=$kV kW=$kW $(suffix)"
            elseif nph == "1"
                bus, phase = split(bus_def, ".")
                load_defs[k] = "$(prefix)Phases=$nph Bus1=$(bus) kV=$kV kW=$kW $(suffix)"
            else
                error("Encountered load with $nph phases.")
            end
        end
    end
    _add_lines_with_header!(lines_load, load_defs, "Imported from Loads.txt, with modified neutral terminal")

    open(loads_file, "w") do f
        write(f, join(lines_load, "\n"))
    end

    return nothing
end


# load the monitor definitions
function build_monitors_txt(case_dir, output_folder)
    monitors_file = "$output_folder/Monitors.txt"
    lines_monitor = Vector{String}()
    monitor_defs = open("$case_dir/Monitors.txt", "r") do f
        readlines(f)
    end
    push!(monitor_defs, "New Monitor.LINE1_1_PQ_vs_Time Line.LINE1 1 Mode=1 ppolar=0")
    push!(monitor_defs, "New Monitor.LINE1_1_VI_vs_Time Line.LINE1 1 Mode=0")
    push!(monitor_defs, "New Monitor.LINE0_PQ_vs_Time Line.LINE0 1 Mode=1 ppolar=0")
    push!(monitor_defs, "New Monitor.LINE0_VI_vs_Time Line.LINE0 1 Mode=0")
    _add_lines_with_header!(lines_monitor, monitor_defs, "Imported from Monitors.txt")

    open(monitors_file, "w") do f
        write(f, join(lines_monitor, "\n"))
    end
    return nothing
end


# build the dss network
function build_case(case_dir, target_lcs, output_folder, n_wire, cct_name, embedded_load, with_transformer)
    build_master_dss(output_folder, cct_name, n_wire, embedded_load, with_transformer)
    if with_transformer
        build_transformers_txt(case_dir, output_folder, n_wire)
    end
    lc_mapping = build_linecode_txt(case_dir, target_lcs, output_folder)
    build_lines_txt(case_dir, output_folder, target_lcs, lc_mapping, n_wire, with_transformer)
    if !embedded_load.is_embedded
        build_loadshapes_txt(case_dir, output_folder)
        build_monitors_txt(case_dir, output_folder)
    end
    build_loads_txt(case_dir, output_folder, n_wire, embedded_load)
end
