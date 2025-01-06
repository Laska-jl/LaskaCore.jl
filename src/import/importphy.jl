# Import output from phy as a PhyOutput struct


"""
    importphy(phydir::String, glxdir::String="", triggerpath::String=""; filter::Union{Pair,Vector{Pair{Symbol,Function}},Nothing}=nothing, includemua::Bool=false, samplerate::Union{Nothing,Real}=nothing)

Import Kilosort output processed in Phy. Spiketimes are sorted.

By default, only "good" clusters as per phy output are included. Setting `includemua=true` will include "mua" clusters as well as unclassified.            

# Arguments

- **Positional arguments**
    - `phydir`: Phy output directory. Has to contain the files "spike_clusters.npy", "spike_times.npy" and "cluster_info.tsv".
    - `glxdir`: SpikeGLX output directory. Has to contain a ".meta" file. May be omitted if the keyword argument `samplerate` is provided. 
    - `triggerpath`: Direct path to a ".csv" or ".bin" file containing a *single* trigger channel as exported from SpikeGLX. Optional argument.
- **Keyword arguments**
    - `filter`: A `Pair` or `Vector{Pair}` allowing clusters to be excluded based on variables in *cluster_info.tsv* (see below for examples).
    - `samplerate`: The samplerate of the spike times. Required if `glxdir` is not specified. If both `glxdir` and `samplerate` is specified the `samplerate` argument will override the samplerate found in the SpikeGLX .meta file.
    - `includemua`: If `true` clusters not marked as good in KiloSort/Phy will be included. `false` by default.


Clusters may be further filtered based on any variable in "cluster_info.tsv". 
This can be done by passing a `Pair{Symbol, Function}` containing the column to filter by and a function to be applied to aforementioned column (see below for examples).
The syntax of the filter itself is exactly the same as when filtering a DataFrame, ie `:somevariable => x -> x > 42` in order to only include clusters where `somevariable` in cluster_info.tsv is greater than 42.
Several filters may be applied by wrapping them in a `Vector`. For further information see the `DataFrames.jl` [documentation on filter](https://dataframes.juliadata.org/stable/lib/functions/#Base.filter).

# Example 
```Julia
### Exclude any cluster with a mean firing rate less than 1:

#                            A `Symbol` identifying the column to apply the filtering function to ----vvv
result = importphy("phyoutput_directory", "glxoutput_directory", "direct_path_to_triggerfile"; filter=:fr => x -> x > 1)
#                                                        An anonymous function returning `true` if x > 1 ----^^^^^^^^^^


### Exclude any cluster with a mean firing rate less than 1 and/or amp less than 30

filtervec = [:fr => x -> x > 1, :amp => x -> x >= 30]

result2 = importphy("phyoutput_directory", "glxoutput_directory", "direct_path_to_triggerfile"; filter=filtervec)

```

# Returns

A [`LaskaCore.PhyOutput`](@ref) struct.

"""
function importphy(phydir::String, glxdir::String="", triggerpath::String=""; filter::Union{Pair,Vector{Pair{Symbol,Function}},Nothing}=nothing, includemua::Bool=false, samplerate::Union{Nothing,Real}=nothing)

    if glxdir == "" && isnothing(samplerate)
        throw(ArgumentError("If no spikeGLX directory with a .meta file is provided the sample rate of spiketimes need to be specified"))
    end

    clusterpath = joinpath(phydir, "spike_clusters.npy")
    amppath = joinpath(phydir, "amplitudes.npy")
    timespath = joinpath(phydir, "spike_times.npy")
    infopath = joinpath(phydir, "cluster_info.tsv")

    clusters, times, info, amps = _importclusterstimesinfoamp(clusterpath, timespath, infopath, amppath)

    idvec = info[!, "cluster_id"] |> vec

    resdict = _sortspiketimes(clusters, times, idvec, amps)

    if !includemua
        isgood(group) = group == "good"
        info = subset(info, :group => ByRow(isgood), skipmissing=true)
    end

    # Filter info df
    if !isnothing(filter)
        if filter isa Vector
            for f in filter
                filter!(f, info)
            end
        else
            filter!(filter, info)
        end
    end

    idvec = deepcopy(info[!, "cluster_id"])

    if triggerpath != ""
        triggers = _importtriggertimes(triggerpath)
    else
        triggers = nothing
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
        metadict = Dict{String,String}(String(i[1]) => String(i[2]) for i in metaraw)
        samprate = parse(Float64, metadict["imSampRate"])
        samprate = isinteger(samprate) ? Int64(samprate) : samprate
        if !isnothing(samplerate) # If sample rate provided as keyword arg, override value found in SpikeGLX.meta file
            @info "Sample rate of spike times set to $samplerate as specified instead of $samprate found in SpikeGLX .meta file"
            samprate = samplerate
        end
    else # If no glxdir provided, set samprate to manually specified value and initialize an empty metadict
        samprate = samplerate
        metadict = nothing
    end

    clustervec = Vector{Cluster{eltype(times),typeof(samprate), eltype(amps)}}(undef, 0)

    __populateclustervec!(clustervec, idvec, info, samprate, resdict)

    # for id in idvec
    #     inf = info[findall(x -> x == id, info[!, "cluster_id"]), :]
    #     push!(clustervec, Cluster(id, inf, SpikeVector(sort!(resdict[id]), samprate)))
    # end

    return PhyOutput(idvec, clustervec, triggers, metadict, info, phydir)
end
