# Returns a dict with keys from 'ids' and values 
# of vectors of spiketimes sorted according to their cluster.
function _sortspiketimes(clusters::Vector{C}, times::Vector{T}, ids::Vector{I}) where {C,T,I}
    d = Dict{C,Vector{T}}(C(c) => Vector{T}(undef, 0) for c in ids)
    _populateresdict!(d, clusters, times)
    return d
end

function _populateresdict!(resdict::Dict{C,Vector{T}}, clusters::Vector{C}, times::Vector{T}) where {C,T}
    length(clusters) != length(times) && throw(AssertionError("'clusters' and 'times' vector not of equal length"))
    @inbounds for i in eachindex(clusters)
        push!(resdict[clusters[i]], times[i])
    end
end

# """
#     _importclusterstimesinfo(clusterpath::String, timespath::String, infopath::String)
#
# Utility function for importing clusters, spiketimes, and info .tsv-file from phy directory.
#
# If spiketimes or cluster ids are a subtype of `Integer` but not an `Int64` they will be converted to `Int64`.
# """
function _importclusterstimesinfo(clusterpath::String, timespath::String, infopath::String)
    clusters = vec(NPZ.npzread(clusterpath))
    times = vec(NPZ.npzread(timespath))
    info = CSV.read(infopath, DataFrame)
    if eltype(clusters) <: Integer && !(eltype(times) isa Int64)
        clusters = convert(Vector{Int64}, clusters)
    end
    if eltype(times) <: Integer && !(eltype(times) isa Int64)
        times = convert(Vector{Int64}, times)
    end
    return clusters, times, info
end

function _importtriggertimes(path::String)
    path == "" && return []
    if path[end-3:end] == ".bin"
        t = importchanint16bin(path)
    elseif path[end-3:end] == ".csv"
        t = importchanint16csv(path)
    else
        throw(ArgumentError(
            """
            Path to triggertimes should be one of:
            - An empty string if no trigger path should be imported
            - A direct path to a .bin or .csv file as exported by SpikeGLX containing a single channel
            """
        ))
    end
    return gettrig(t)
end
