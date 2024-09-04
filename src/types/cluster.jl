####################################################
#
# Cluster is meant to hold 1 cluster as well as all
# info such as depth, fr and computed values such as MAD anc cv2.
#
####################################################


"""
    abstract type AbstractCluster{T,U} end

Parent type to concrete types representing single clusters.

"""
abstract type AbstractCluster{T,U} end



"""
    struct Cluster{T,U} <: AbstractCluster{T,U}
        id::Int64
        info::DataFrame
        spiketimes::SpikeVector{T,U}
    end

Struct for holding a single Cluster.

Direct field access is **not** recommended. Basic interface functions include:

- [`LaskaCore.id`](@ref) -- Returns the Cluster id.
- [`LaskaCore.nspikes`](@ref) -- Returns the number of spikes in the cluster (Based off length of the `spiketimes` field).
- [`LaskaCore.info`](@ref) -- Returns the info of the `Cluster` from "cluster_info.tsv" as a `DataFrame` row.
- [`LaskaCore.spiketimes`](@ref) -- Returns a Vector containing all spiketimes of the `Cluster`.


"""
struct Cluster{T,U} <: AbstractCluster{T,U}
    id::Int64
    info::DataFrame
    spiketimes::SpikeVector{T,U}
end


"""
    id(cluster::T) where {T<:AbstractCluster}

Returns the id of `cluster`
"""
function id(cluster::T) where {T<:AbstractCluster}
    return cluster.id
end

"""
    nspikes(cluster::T) where {T<:AbstractCluster}

Returns the number of spikes of `cluster`.

"""
function nspikes(cluster::Cluster)
    return length(cluster.spiketimes)
end

"""
    info(cluster::T) where {T<:AbstractCluster}
    info(cluster::T, var::String) where {T<:AbstractCluster}

Returns info (as a DataFrame) about `cluster`. A string may be supplied to return a specific entry.
"""
function info(cluster::T) where {T<:AbstractCluster}
    return cluster.info
end

function info(cluster::T, var::String) where {T<:AbstractCluster}
    return cluster.info[1, var]
end



"""
    spiketimes(cluster::Cluster::T) where {T<:AbstractCluster}

Returns the spiketimes of `cluster`.
"""
function spiketimes(cluster::T) where {T<:AbstractCluster}
    return cluster.spiketimes
end

"""
    samplerate(cluster::AbstractCluster)

Returns the current samplerate (in Hz) of `cluster`.
"""
function samplerate(cluster::AbstractCluster)
    return spiketimes(cluster) |> samplerate
end




"""
    struct RelativeCluster{T,U} <: AbstractCluster{T,U}
        id::Int64
        info::SubDataFrame
        spiketimes::RelativeSpikeVector{T,U}
    end


Struct for holding a cluster and its spiketimes relative to triggers.       
Similar to `Cluster{T}` except that the field `spiketimes` is a `Vector{Vector{T}}` where each vector represents trigger #n.

Direct field access is **not** recommended. Basic interface functions include:

- [`LaskaCore.id`](@ref) -- Returns the Cluster id.
- [`LaskaCore.nspikes`](@ref) -- Returns the number of spikes in each spiketrain.
- [`LaskaCore.info`](@ref) -- Returns the info of the `Cluster` from "cluster_info.tsv" as a `DataFrame` row.
- [`LaskaCore.spiketimes`](@ref) -- Returns a Vector containing all spiketimes of the `Cluster`.
- [`LaskaCore.samplerate`](@ref) -- Returns the samplerate of the spiketimes belonging to the cluster.

"""
struct RelativeCluster{T,U} <: AbstractCluster{T,U}
    id::Int64
    info::DataFrame
    spiketimes::RelativeSpikeVector{T,U}
end

"""
    nspikes(cluster::RelativeCluster)

Returns the total number of spikes around each trigger in `cluster`.
"""
function nspikes(cluster::RelativeCluster)
    return length.(spiketimes(cluster))
end

function Base.show(io::IO, obj::Cluster)
    println("$(typeof(obj)) with id=$(id(obj)) and $(nspikes(obj)) spikes")
end

function Base.show(io::IO, obj::RelativeCluster)
    println("$(typeof(obj)) with id=$(id(obj)) and $(nspikes(obj)) spikes across $(length(spiketimes(obj))) trigger events")
end
