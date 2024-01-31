####################################################
#
# Cluster is meant to hold 1 cluster as well as all
# info such as depth, fr and computed values such as MAD anc cv2.
#
####################################################


"""
    abstract type AbstractCluster{T} end

Parent type to concrete types representing single clusters.

"""
abstract type AbstractCluster{T} end



"""
    struct Cluster{T} <: AbstractCluster{T}
        id::Int64
        info::SubDataFrame
        spiketimes::Vector{T}
    end

Struct for holding a single Cluster.

Direct field access is **not** recommended. Basic interface functions include:

- [`LaskaCore.id`](@ref) -- Returns the Cluster id.
- [`LaskaCore.nspikes`](@ref) -- Returns the number of spikes in the cluster (Based off length of the `spiketimes` field).
- [`LaskaCore.info`](@ref) -- Returns the info of the `Cluster` from "cluster_info.tsv" as a `DataFrame` row. This is a *view* of the DataFrame in the parent `Experiment`.
- [`LaskaCore.spiketimes`](@ref) -- Returns a Vector containing all spiketimes of the `Cluster`.


"""
struct Cluster{T} <: AbstractCluster{T}
    id::Int64
    info::DataFrame
    spiketimes::Vector{T}
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

# sample frequency/time conversion


"""
    timetosamplerate(cluster::T, time::U) where {T<:AbstractCluster, U<:Quantity{<:Number, Unitful.𝐓}}

Convert a `time` of some unit from `Unitful.jl` into the samplerate of the `cluster`. For example, converting 10ms to a samplerate of 30 000Hz will yield a value of 300.0.

# Examples

```julia
using LaskaCore
using Unitful

c = getcluster(exp, 33) # Get cluster '33' from an AbstractExperiment

timetosamplerate(c, 10u"ms") # Convert 10ms into the samplerate of the cluster
```
"""
function timetosamplerate(cluster::T, time::U) where {T<:AbstractCluster, U<:Quantity{<:Number, Unitful.𝐓}}
	samp::Float64 = @views info(cluster, "samprate")
    samp * ustrip(u"s", time)
end

function timetosamplerate(cluster::T, time::U) where {T<:AbstractCluster, U<:AbstractRange{<:Quantity{<:Number, Unitful.𝐓}}}
	samp::Float64 = @views info(cluster, "samprate")
    (samp * ustrip(u"s", time[begin])):(samp * ustrip(u"s", time[end]))

end
# Indexing

function Base.getindex(cluster::T, I::U) where {T<:AbstractCluster, U<:Integer}
    1 <= I <= length(spiketimes(cluster)) || throw(BoundsError(cluster, I))
	@inline spiketimes(cluster)[I]
end

function Base.getindex(cluster::T, I::U) where {T<:AbstractCluster, U<:Number}
    cluster[convert(Int, I)]
end

function Base.getindex(cluster::T, I::Vector{U}) where {T<:AbstractCluster, U<:Number}
	[cluster[i] for i in I]
end


Base.firstindex(cluster::T) where {T<:AbstractCluster} = 1
Base.lastindex(cluster::T) where {T<:AbstractCluster} = length(spiketimes(cluster))

# Time unit based indexing

"""
    Base.getindex(cluster::T, I::U) where {T<:AbstractCluster, U<:AbstractRange{<:Quantity{<:Number, Unitful.𝐓}}}

Allows indexing with time ranges from `Unitful.jl`.
"""
function Base.getindex(cluster::T, I::U) where {T<:AbstractCluster, U<:AbstractRange{<:Quantity{<:Number, Unitful.𝐓}}}
    rang = timetosamplerate(cluster, I)
    filter(x -> rang[begin] <= x <= rang[end], spiketimes(cluster))
end

"""
    struct RelativeCluster{T} <: AbstractCluster{T}
        id::Int64
        info::SubDataFrame
        spiketimes::Vector{Vector{T}}
    end


Struct for holding a cluster and its spiketimes relative to triggers.       
Similar to `Cluster{T}` except that the field `spiketimes` is a `Vector{Vector{T}}` where each vector represents trigger #n.

Direct field access is **not** recommended. Basic interface functions include:

- [`LaskaCore.id`](@ref) -- Returns the Cluster id.
- [`LaskaCore.nspikes`](@ref) -- Returns the number of spikes in the cluster (Based off length of the `spiketimes` field).
- [`LaskaCore.info`](@ref) -- Returns the info of the `Cluster` from "cluster_info.tsv" as a `SubDataFrame`. This is a *view* of the `info` `DataFrame` from the parent `Experiment`.
- [`LaskaCore.spiketimes`](@ref) -- Returns a Vector containing all spiketimes of the `Cluster`.

"""
struct RelativeCluster{T} <: AbstractCluster{T}
    id::Int64
    info::DataFrame
    spiketimes::Vector{Vector{T}}
end
