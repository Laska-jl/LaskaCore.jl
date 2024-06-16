#-------------------------------------------------------------
# Container(s) for all Clusters that are part of an experiment
#-------------------------------------------------------------


"""
    abstract type AbstractExperiment{T} end

Parent type to concrete types representing entire experiments containing their specifications, metadata and clusters.
"""
abstract type AbstractExperiment{T,U,M<:Union{Nothing,Dict{String,String}},S<:Union{Nothing,Vector{<:Integer}}} end

"""
    struct PhyOutput{T} <: AbstractExperiment{T}
        clusterids::Vector{Int64}
        clusters::Vector{Cluster{T}}
        trigtimes::Vector{T}
        meta::Dict{SubString{String},SubString{String}}
        info::DataFrame
    end

Struct for holding Kilosort output preprocessed in Phy. Should be instantiated using the outer constructor [`LaskaCore.importphy`](@ref).

Direct field access is **not** recommended. Basic interface functions include:

- [`LaskaCore.clusterids`](@ref) -- Returns all cluster ID:s as a Vector.
- [`LaskaCore.getcluster`](@ref) -- Returns a specific [`LaskaCore.Cluster`](@ref).
- [`LaskaCore.clustervector`](@ref) -- Returns all [`LaskaCore.Cluster`](@ref):s.
- [`LaskaCore.getmeta`](@ref) -- Returns the spikeGLX meta as a dict or a specific entry.
- [`LaskaCore.triggertimes`](@ref) -- Returns the trigger event times.
- [`LaskaCore.ntrigs`](@ref) -- Returns the length of the trigger event time Vector.

"""
struct PhyOutput{T,U,M,S} <: AbstractExperiment{T,U,M,S}
    clusterids::Vector{Int64}
    clusters::Vector{Cluster{T,U}}
    trigtimes::S
    meta::M
    info::DataFrame
end


"""
    struct RelativeSpikes{T} <: AbstractExperiment{T}
        clusterids::Vector{Int64}
        clusters::Vector{RelativeCluster{T}}
        trigtimes::Vector{T}
        meta::Dict{SubString{String},SubString{String}}
        info::DataFrame
        stimtrain::Dict{String,T}
        specs::Dict{String,T}
    end

Similar to [`LaskaCore.PhyOutput`](@ref). However, instead of [`LaskaCore.Cluster`](@ref)s, [`LaskaCore.RelativeCluster`](@ref)s are used.
In these, spiketimes relative to trigger event(s) instead of absolute spiketimes are used.              
Additionally contains the fields:

- `stimtrain::Dict{String,T}` -- Dict in the format `name/id` => `time`. `time` should be relative to trigger event. Specified by the user on creation of struct using [`LaskaCore.relativespikes`](@ref).
- `specs::Dict{String,T}` -- Dict containing the time before/after (`back`/`forward`) trigger events that spikes are included; as well as number of trigger events (`ntrig`).

"""
struct RelativeSpikes{T,U,M,S,V} <: AbstractExperiment{T,U,M,S}
    clusterids::Vector{Int64}
    clusters::Vector{RelativeCluster{T,U}}
    trigtimes::S
    meta::M
    info::DataFrame
    stimtrain::Dict{String,T}
    specs::Dict{String,V}
end


"""
    getcluster(experiment::T, cluster::Int) where {T<:AbstractExperiment}


Returns a `cluster` from `experiment`.

"""
function getcluster(experiment::T, cluster) where {T<:AbstractExperiment}
    return experiment.clusters[findfirst(x -> x == cluster, clusterids(experiment))]
end


"""

    ntrigs(experiment::T) where {T<:AbstractExperiment}

Returns the number of trigger events in `experiment`.
"""
function ntrigs(experiment::T) where {T<:AbstractExperiment}
    return length(triggertimes(experiment))
end

"""
    clusterids(experiment::T) where {T<:AbstractExperiment}

Returns a Vector of all cluster id:s present in experiment.

"""
function clusterids(experiment::T) where {T<:AbstractExperiment}
    return experiment.clusterids
end

"""
    triggertimes(experiment::T) where {T<:AbstractExperiment}

Returns the timestamps of trigger events in `experiment`.

"""
function triggertimes(experiment::T) where {T<:AbstractExperiment}
    return experiment.trigtimes
end

"""
    clustervector(experiment::T) where {T<:AbstractExperiment}

Returns a `Vector{T}` where T<:AbstractCluster containing all clusters in `experiment`.

"""
@views function clustervector(experiment::T) where {T<:AbstractExperiment}
    return experiment.clusters
end



"""
    spiketimes(experiment::PhyOutput)

Get all spiketimes in `experiment`. Spiketimes are **not sorted** by time.
"""
function spiketimes(experiment::PhyOutput)
    vcat(spiketimes.(clustervector(experiment))...)
end

"""
    getmeta(experiment::T, entry::String) where {T<:AbstractExperiment}
    getmeta(experiment::T) where {T<:AbstractExperiment}

Returns experiment meta info from spikeGLX. If an `entry` string is not supplied all entries are returned.

"""
function getmeta end

function getmeta(experiment::T, entry::String) where {T<:AbstractExperiment}
    return experiment.meta[entry]
end

function getmeta(experiment::T) where {T<:AbstractExperiment}
    return experiment.meta
end

"""
    info(experiment::T) where {T<:AbstractExperiment}
    info(experiment::T, var) where {T<:AbstractExperiment}

Returns the `cluster_info.tsv` attached to the `experiment` in the form of a `DataFrame`.
If `var` is provided returns a `Vector` of the matching column.
"""
function info end

function info(experiment::T) where {T<:AbstractExperiment}
    return experiment.info
end

function info(experiment::T, var) where {T<:AbstractExperiment}
    return experiment.info[!, var]
end

# RelativeSpikes-specific functions
"""
    relativespecs(rel::RelativeSpikes{T}) where {T<:Real}
    relativespecs(rel::RelativeSpikes{T}, spec::String) where {T<:Real}


Returns a Dict containing the 'specs' of a `RelativeSpikes` struct.                         
Includes the `back` and `forward` variables used as well as the number of trigger events (`ntrigs`).
"""
function relativespecs end

function relativespecs(rel::RelativeSpikes{T}) where {T<:Real}
    return rel.specs
end


function relativespecs(rel::RelativeSpikes{T}, spec::String) where {T<:Real}
    return rel.specs[spec]
end

"""
    stimtimes(experiment::RelativeSpikes)

Returns a dict containing the stimtrain of a `RelativeSpikes` struct in the form of `label => time`.
"""
function stimtimes(experiment::RelativeSpikes)
    return experiment.stimtrain
end

"""
    function spiketimes(experiment::RelativeSpikes{T}) where T

Get all spiketimes of `experiment`.
"""
function spiketimes(experiment::RelativeSpikes{T}) where {T}
    Nt = ntrigs(experiment)
    outvec = [T[] for _ in 1:Nt]
    for c in spiketimes.(clustervector(experiment))
        for (n, v) in enumerate(c)
            push!(outvec[n], v...)
        end
    end
    return outvec
end

"""
    nclusters(exp::AbstractExperiment)

Returns the number of clusters in `exp`.
"""
nclusters(exp::AbstractExperiment) = length(clusterids(exp))

# Show and such

function Base.show(io::IO, obj::AbstractExperiment)
    println("$(typeof(obj)) containing $(nclusters(obj)) $(typeof(clustervector(obj)[1])):\n$(clusterids(obj))")
end


