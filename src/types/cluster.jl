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
        info::SubDataFrame
spiketimes::SpikeVector{T,U}
    end

Struct for holding a single Cluster.

Direct field access is **not** recommended. Basic interface functions include:

- [`LaskaCore.id`](@ref) -- Returns the Cluster id.
- [`LaskaCore.nspikes`](@ref) -- Returns the number of spikes in the cluster (Based off length of the `spiketimes` field).
- [`LaskaCore.info`](@ref) -- Returns the info of the `Cluster` from "cluster_info.tsv" as a `DataFrame` row. This is a *view* of the DataFrame in the parent `Experiment`.
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

Returns the number of spikes (length of the `spiketimes` field) of `Cluster`.

"""
function nspikes(cluster::T) where {T<:AbstractCluster}
    return length(cluster.spiketimes)
end

"""
    info(cluster::T) where {T<:AbstractCluster}

Returns info (as dict) about `cluster`. A string may be supplied to return a specific entry (as Float64).
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

function samplerate(cluster::AbstractCluster)
    return spiketimes(cluster) |> samplerate
end

# Indexing

function Base.getindex(cluster::T, I::U) where {T<:AbstractCluster,U<:Integer}
    1 <= I <= length(spiketimes(cluster)) || throw(BoundsError(cluster, I))
    @inline spiketimes(cluster)[I]
end

function Base.getindex(cluster::T, I::U) where {T<:AbstractCluster,U<:Number}
    cluster[convert(Int, I)]
end

function Base.getindex(cluster::T, I::Vector{U}) where {T<:AbstractCluster,U<:Number}
    [cluster[i] for i in I]
end


Base.firstindex(cluster::T) where {T<:AbstractCluster} = 1
Base.lastindex(cluster::T) where {T<:AbstractCluster} = length(spiketimes(cluster))

# Time unit based indexing

"""
    Base.getindex(cluster::T, I::U) where {T<:AbstractCluster, U<:AbstractRange{<:Quantity{<:Number, Unitful.𝐓}}}

Allows indexing with time ranges from `Unitful.jl`.
"""
function Base.getindex(cluster::T, I::U) where {T<:AbstractCluster,U<:AbstractRange{<:Quantity{<:Number,Unitful.𝐓}}}
    rang = timetosamplerate(cluster, I)
    filter(x -> rang[begin] <= x <= rang[end], spiketimes(cluster))
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
- [`LaskaCore.nspikes`](@ref) -- Returns the number of spikes in the cluster (Based off length of the `spiketimes` field).
- [`LaskaCore.info`](@ref) -- Returns the info of the `Cluster` from "cluster_info.tsv" as a `SubDataFrame`. This is a *view* of the `info` `DataFrame` from the parent `Experiment`.
- [`LaskaCore.spiketimes`](@ref) -- Returns a Vector containing all spiketimes of the `Cluster`.

"""
struct RelativeCluster{T,U} <: AbstractCluster{T,U}
    id::Int64
    info::DataFrame
    spiketimes::RelativeSpikeVector{T,U}
end
