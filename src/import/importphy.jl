# Import output from phy as a PhyOutput struct

#  TODO: Add separate docstrings
"""

    importphy(phydir::String, glxdir::String, triggerpath::String; includemua::Bool=false)
    importphy(phydir::String, filters::Tuple{Symbol,Function}, glxdir::String="", triggerpath::String=""; includemua::Bool=false)
    importphy(phydir::String, filters::Vector{Tuple{Symbol,Function}}, glxdir::String="", triggerpath::String=""; includemua::Bool=false)

Import Kilosort output processed in Phy. Spiketimes are sorted.

By default, only "good" clusters as per phy output are included. Setting `includemua=true` will include "mua" clusters as well as unclassified.            
Clusters may be further filtered based on any variable in "cluster_info.tsv". This is done by including a Tuple
with the column to be filtered as a symbol and a filtering function. Several Tuples containing such filters may be included by wrapping them in a `Vector`.         



# Example 
```Julia
# Removes any cluster with a mean firing rate less than 1:

# A function returning true if x > 1
filterfr = x -> x > 1
filteramp = x -> x >= 100

# A tuple with the above function (filterfunc) and the column to which it should be
# applied (:fr).
frtuple = (:fr, filterfr)

# A second tuple with a function returning `true` if x >= 100 and the column (:amp) to which
# it should be applied
amptuple = (:amp, filteramp)

# A `Vector` with the 2 `Tuple`s above.
filtervec = [frtuple, amptuple]

result = importphy("phyoutput_directory", filtervec, "glxoutput_directory", "direct_path_to_triggerfile")
```

"""
function importphy(phydir::String, glxdir::String="", triggerpath::String=""; includemua::Bool=false)
    if Sys.iswindows()
        clusters::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(phydir * "\\spike_clusters.npy"))
        times::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(phydir * "\\spike_times.npy")[:, 1])
        info::DataFrames.DataFrame = CSV.read(phydir * "\\cluster_info.tsv", DataFrame)
    else
        clusters = convert(Vector{Int64}, NPZ.npzread(phydir * "/spike_clusters.npy"))
        times = convert(Vector{Int64}, NPZ.npzread(phydir * "/spike_times.npy")[:, 1])
        info = CSV.read(phydir * "/cluster_info.tsv", DataFrame)
    end

    idvec::Vector{Int64} = info[!, "cluster_id"]

    clustervec::Vector{Cluster{Int64}} = Vector(undef, 0)

    resdict::Dict{Int64,Vector{Int64}} = Dict(c => Vector{Int64}(undef, 0) for c in idvec)

    for _ in eachindex(clusters)
        cluster::Int64 = pop!(clusters)
        time::Int64 = pop!(times)
        push!(resdict[cluster], time)
    end
    if !includemua
        isgood(group) = group == "good"
        info = subset(info, :group => ByRow(isgood), skipmissing=true)
    end
    idvec = deepcopy(info[!, "cluster_id"])


    if triggerpath != ""
        if triggerpath[end-3:end] == ".bin"
            t = importchanint16bin(triggerpath)
        elseif triggerpath[end-3:end] == ".csv"
            t = importchanint16csv(triggerpath)
        else
            error("Triggerfile (triggerpath) must end in '.csv' or 'bin'!")
        end
        triggers = gettrig(t)
    else
        triggers = Vector{Int64}(undef, 0)
    end

    if glxdir != ""
        glxfiles = readdir(glxdir, join=true)
        binlist = [f for f in glxfiles if f[Base.length(f)-6:Base.length(f)] == ".ap.bin"]
        if length(binlist) > 0
            binfile::String = binlist[1]
        else
            binfile = ""
        end
        metafile::String = [f for f in glxfiles if f[Base.length(f)-7:Base.length(f)] == ".ap.meta"][1]

        # Read metadata
        tmp = open(metafile, "r")
        metaraw = readlines(tmp)
        close(tmp)
        metaraw = split.(metaraw, "=")
        metadict = Dict{SubString{String},SubString{String}}(i[1] => i[2] for i in metaraw)
        samprate = parse(Float64, metadict["imSampRate"])
        samprate = isinteger(samprate) ? Int64(samprate) : samprate
        info.samprate = [samprate for _ in 1:size(info, 1)]
    else
        metadict = Dict{SubString{String},SubString{String}}()
    end


    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, SpikeVector(sort!(resdict[id]), samprate)))
    end

    return PhyOutput(idvec, clustervec, triggers, metadict, info)
end


# Version with filters
function importphy(phydir::String, filters::Tuple{Symbol,Function}, glxdir::String="", triggerpath::String=""; includemua::Bool=false)
    if Sys.iswindows()
        clusters::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(phydir * "\\spike_clusters.npy"))
        times::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(phydir * "\\spike_times.npy")[:, 1])
        info::DataFrames.DataFrame = CSV.read(phydir * "\\cluster_info.tsv", DataFrame)
    else
        clusters = convert(Vector{Int64}, NPZ.npzread(phydir * "/spike_clusters.npy"))
        times = convert(Vector{Int64}, NPZ.npzread(phydir * "/spike_times.npy")[:, 1])
        info = CSV.read(phydir * "/cluster_info.tsv", DataFrame)
    end

    idvec::Vector{Int64} = info[!, "cluster_id"]

    clustervec::Vector{Cluster{Int64}} = Vector(undef, 0)

    resdict::Dict{Int64,Vector{Int64}} = Dict(c => Vector{Int64}(undef, 0) for c in idvec)

    for _ in eachindex(clusters)
        cluster::Int64 = pop!(clusters)
        time::Int64 = pop!(times)
        push!(resdict[cluster], time)

    end

    if !includemua
        isgood(group) = group == "good"
        info = subset(info, :group => ByRow(isgood), skipmissing=true)
    end
    filter!(filters[1] => filters[2], info)
    idvec = deepcopy(info[!, "cluster_id"])

    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, sort!(resdict[id])))
    end

    if triggerpath != ""
        t = importchanint16(triggerpath)
        triggers = gettrig(t)
    else
        triggers = Vector{Int64}(undef, 0)
    end

    if glxdir != ""
        glxfiles = readdir(glxdir, join=true)
        binlist = [f for f in glxfiles if f[Base.length(f)-6:Base.length(f)] == ".ap.bin"]
        if length(binlist) > 0
            binfile::String = binlist[1]
        else
            binfile = ""
        end
        metafile::String = [f for f in glxfiles if f[Base.length(f)-7:Base.length(f)] == ".ap.meta"][1]

        # Read metadata
        tmp = open(metafile, "r")
        metaraw = readlines(tmp)
        close(tmp)
        metaraw = split.(metaraw, "=")
        metadict = Dict{SubString{String},SubString{String}}(i[1] => i[2] for i in metaraw)
        samprate = parse(Float64, metadict["imSampRate"])
        samprate = isinteger(samprate) ? Int64(samprate) : samprate
        #TODO: Remove this since sample rate is now contained in the spike vectors themselves
        info.samprate = [samprate for _ in 1:size(info, 1)]
    else
        metadict = Dict{SubString{String},SubString{String}}()
    end


    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, sort!(resdict[id]), samprate))
    end

    return PhyOutput(idvec, clustervec, triggers, metadict, info)
end


# Import with several filters in the form of a vector
function importphy(phydir::String, filters::Vector{Tuple{Symbol,Function}}, glxdir::String="", triggerpath::String=""; includemua::Bool=false)
    if Sys.iswindows()
        clusters::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(phydir * "\\spike_clusters.npy"))
        times::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(phydir * "\\spike_times.npy")[:, 1])
        info::DataFrames.DataFrame = CSV.read(phydir * "\\cluster_info.tsv", DataFrame)
    else
        clusters = convert(Vector{Int64}, NPZ.npzread(phydir * "/spike_clusters.npy"))
        times = convert(Vector{Int64}, NPZ.npzread(phydir * "/spike_times.npy")[:, 1])
        info = CSV.read(phydir * "/cluster_info.tsv", DataFrame)
    end

    idvec::Vector{Int64} = info[!, "cluster_id"]

    clustervec::Vector{Cluster{Int64}} = Vector(undef, 0)

    resdict::Dict{Int64,Vector{Int64}} = Dict(c => Vector{Int64}(undef, 0) for c in idvec)

    for _ in eachindex(clusters)
        cluster::Int64 = pop!(clusters)
        time::Int64 = pop!(times)
        push!(resdict[cluster], time)

    end

    if !includemua
        isgood(group) = group == "good"
        info = subset(info, :group => ByRow(isgood), skipmissing=true)
    end

    for f in filters
        filter!(f[1] => f[2], info)
    end
    idvec = deepcopy(info[!, "cluster_id"])


    if triggerpath != ""
        t = importchanint16(triggerpath)
        triggers = gettrig(t)
    else
        triggers = Vector{Int64}(undef, 0)
    end

    if glxdir != ""
        glxfiles = readdir(glxdir, join=true)
        binlist = [f for f in glxfiles if f[Base.length(f)-6:Base.length(f)] == ".ap.bin"]
        if length(binlist) > 0
            binfile::String = binlist[1]
        else
            binfile = ""
        end
        metafile::String = [f for f in glxfiles if f[Base.length(f)-7:Base.length(f)] == ".ap.meta"][1]

        # Read metadata
        tmp = open(metafile, "r")
        metaraw = readlines(tmp)
        close(tmp)
        metaraw = split.(metaraw, "=")
        metadict = Dict{SubString{String},SubString{String}}(i[1] => i[2] for i in metaraw)
        samprate = parse(Float64, metadict["imSampRate"])
        samprate = isinteger(samprate) ? Int64(samprate) : samprate
        info.samprate = [samprate for _ in 1:size(info, 1)]
    else
        metadict = Dict{SubString{String},SubString{String}}()
    end

    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, sort!(resdict[id]), samprate))
    end

    return PhyOutput(idvec, clustervec, triggers, metadict, info)
end
