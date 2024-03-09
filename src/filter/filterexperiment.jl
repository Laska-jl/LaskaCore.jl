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
    exclude = Vector{Bool}(undef, nclusters(exp))
    filter!(var => f, exp.info)
    idset = Set(exp.info[:, :cluster_id])
    if isempty(idset)
        @warn "This filter exludes all clusters."
    end
    for (i, c) in enumerate(clusterids(exp))
        exclude[i] = !(c in idset)
    end
    deleteat!(exp.clusters, exclude)
    deleteat!(exp.clusterids, exclude)
    exp
end

"""
    Base.filter((var, f)::Pair, exp::AbstractExperiment{T}) where {T}

Filter the `Cluster`s of `exp` based off a column in its info `DataFrame`. The syntax is the same as when filtering a "normal" DataFrame.

# Example

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
function Base.filter((var, f)::Pair, exp::AbstractExperiment{T}) where {T}
    out = deepcopy(exp)
    filter!(var => f, out)
    return out
end
