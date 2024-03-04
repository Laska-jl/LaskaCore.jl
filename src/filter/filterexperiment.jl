# Functions for filtering experiments

"""
    Base.filter!((var, f)::Pair, exp::AbstractExperiment{T}) where {T}

Filter the `Cluster`s of `exp` based off a column in its info `DataFrame`. The syntax is the same as when filtering a "normal" DataFrame.
`exp` is filtered in-place.

# Example

In order to remove all clusters of `exp` with a depth less than 1000:

```julia
filter!(:depth => x -> x >= 1000.0, exp)
```
"""
function Base.filter!((var, f)::Pair, exp::AbstractExperiment{T}) where {T}
    tmpinfo = filter(var => f, exp.info)
    idset = Set(tmpinfo.cluster_id)
    if isempty(idset)
        throw(ErrorException("No clusters remaining in experiment."))
    end
    idsout = Vector{Int}(undef, size(info(exp), 1))
    clustersout = Vector{Cluster{T}}(undef, nclusters(exp))
    n = 0
    for (id, cluster) in zip(clusterids(exp), clustervector(exp))
        if id in idset
            n += 1
            idsout[n] = id
            clustersout[n] = cluster
        end
    end
    exp.clusterids = idsout[begin:n]
    exp.clusters = clustersout[begin:n]
    exp.info = tmpinfo
end

"""
    Base.filter((var, f)::Pair, exp::PhyOutput{T}) where {T<:Number}

Filter [`Cluster`](@ref)s in a [`PhyOutput`](@ref) based on `fun` which should act on a column of the `info` `DataFrame`.

# Examples

This example will filter a [`PhyOutput`](@ref) containing only [`Cluster`](@ref)s with a contamination percentage less than 5.
```julia
res = importphy(
    PATH_TO_PHYOUTPUT,
    PATH_TO_GLXMETA,
    PATH_TO_TRIGGERCHANNEL
)

filtered = filter(:ContamPct => x -> x <= 5, res)
```

"""
function Base.filter((var, f)::Pair, exp::PhyOutput{T}) where {T<:Number}
    infoout = deepcopy(info(exp))
    filter!(fun, infoout)
    idset = Set(infoout.cluster_id)
    idsout = Vector{Int}(undef, size(info(exp), 1))
    clustersout = Vector{Cluster{T}}(undef, size(info(exp), 1))
    n = 0
    for (id, cluster) in zip(clusterids(exp), clustervector(exp))
        if id in idset
            n += 1
            idsout[n] = id
            clustersout[n] = cluster
        end
    end
    PhyOutput(
        idsout[begin:n],
        clustersout[begin:n],
        triggertimes(exp),
        getmeta(exp),
        infoout
    )
end
