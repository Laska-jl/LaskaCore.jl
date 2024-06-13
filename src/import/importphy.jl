# Import output from phy as a PhyOutput struct

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
function importphy end

function importphy(phydir::String, glxdir::String, triggerpath::String, includemua::Bool=false)
    clusters::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_clusters.npy")))
    times::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_times.npy"))[:, 1])
    info::DataFrames.DataFrame = CSV.read(joinpath(phydir, "cluster_info.tsv"), DataFrame)

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


    if triggerpath[end-3:end] == ".bin"
        t = importchanint16bin(triggerpath)
    elseif triggerpath[end-3:end] == ".csv"
        t = importchanint16csv(triggerpath)
    else
        error("Triggerfile (triggerpath) must be A csv or bin file")
    end
    triggers = gettrig(t)

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
    metadict = Dict{String,String}(String(i[1]) => String(i[2]) for i in metaraw)
    samprate = parse(Float64, metadict["imSampRate"])
    samprate = isinteger(samprate) ? Int64(samprate) : samprate
    info.samprate = [samprate for _ in 1:size(info, 1)]


    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, SpikeVector(sort!(resdict[id]), samprate)))
    end

    return PhyOutput(idvec, clustervec, triggers, metadict, info)
end


# Version with filters
function importphy(phydir::String, glxdir::String, triggerpath::String, filter::Pair, includemua::Bool=false)
    clusters::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_clusters.npy")))
    times::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_times.npy"))[:, 1])
    info::DataFrames.DataFrame = CSV.read(joinpath(phydir, "cluster_info.tsv"), DataFrame)

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

    # Filter info dataframe
    filter!(filter, info)
    idvec = deepcopy(info[!, "cluster_id"])


    t = importchanint16(triggerpath)
    triggers = gettrig(t)

    glxfiles = readdir(glxdir, join=true)
    binlist = [f for f in glxfiles if f[Base.length(f)-6:Base.length(f)] == ".ap.bin"]
    if length(binlist) > 0
        binfile::String = binlist[1]
    else
        binfile = ""
        @warn "Raw spikeGLX file not found (could not find a file ending in '.ap.bin' in specified spikeGLX directory)"
    end

    metafilevec = [f for f in glxfiles if f[Base.length(f)-7:Base.length(f)] == ".ap.meta"]
    if iszero(length(metafilevec))
        throw(ErrorException("Could not find spikeGLX meta file not found (file ending in '.ap.meta' in specified spikeGLX directory)\nIf this file is not available please specify sample rate when calling 'importphy'"))
    elseif length(metafilevec) > 1
        throw(ErrorException("Several files ending in '.ap.meta' found in specified spikeGLX directory"))
    end
    metafile = metafilevec[1]

    # Read metadata
    tmp = open(metafile, "r")
    metaraw = readlines(tmp)
    close(tmp)
    metaraw = split.(metaraw, "=")
    metadict = Dict{String,String}(String(i[1]) => String(i[2]) for i in metaraw)
    samprate = parse(Float64, metadict["imSampRate"])
    samprate = isinteger(samprate) ? Int64(samprate) : samprate
    #TODO: Remove this since sample rate is now contained in the spike vectors themselves
    info.samprate = [samprate for _ in 1:size(info, 1)]

    # Push clusters in idvec which contains filtered clusters
    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, SpikeVector(sort!(resdict[id]), samprate)))
    end

    return PhyOutput(idvec, clustervec, triggers, metadict, info)
end


# Import with several filters in the form of a vector
function importphy(phydir::String, glxdir::String, triggerpath::String, filters::Vector{Pair{Symbol,Function}}; includemua::Bool=false)
    clusters::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_clusters.npy")))
    times::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_times.npy"))[:, 1])
    info::DataFrames.DataFrame = CSV.read(joinpath(phydir, "cluster_info.tsv"), DataFrame)

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
    metadict = Dict{String,String}(String(i[1]) => String(i[2]) for i in metaraw)
        samprate = parse(Float64, metadict["imSampRate"])
        samprate = isinteger(samprate) ? Int64(samprate) : samprate
        info.samprate = [samprate for _ in 1:size(info, 1)]

    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, SpikeVector(sort!(resdict[id]), samprate)))
    end

    return PhyOutput(idvec, clustervec, triggers, metadict, info)
end

# Import without spikeGLX directory

function importphy(phydir::String; triggerpath::String, samplerate::Real, includemua::Bool=false)
    clusters::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_clusters.npy")))
    times::Vector{Int64} = convert(Vector{Int64}, NPZ.npzread(joinpath(phydir, "spike_times.npy"))[:, 1])
    info::DataFrames.DataFrame = CSV.read(joinpath(phydir, "cluster_info.tsv"), DataFrame)

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

    t = importchanint16(triggerpath)
    triggers = gettrig(t)

    for id in idvec
        inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
        push!(clustervec, Cluster(id, inf, SpikeVector(sort!(resdict[id]), samplerate)))
    end

    metadict = Dict{String, String}()

    return PhyOutput(idvec, clustervec, triggers, metadict, info)
    
end
