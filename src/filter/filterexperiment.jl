# Functions for filtering experiments

"""
    Base.filter!(f, exp::AbstractExperiment{T}) where {T}

Filter the `Cluster`s of `exp` based off a column in its info `DataFrame`. The syntax is the same as when filtering a "normal" DataFrame.
`exp` is filtered in-place.

# Example

In order to remove all clusters of `exp` with a depth less than 1000:

```julia
filter!(:depth => x -> x >= 1000.0, exp)
```

For more in-depth information see the documentation for filtering of DataFrames: https://dataframes.juliadata.org/stable/lib/functions/#Base.filter
"""
function Base.filter!(f, exp::AbstractExperiment{T}) where {T}
    exclude = Vector{Bool}(undef, nclusters(exp))
    filter!(f, exp.info)
    idset = Set(exp.info[:, :cluster_id])
    if isempty(idset)
        @warn "This filter exludes all clusters."
    end
    for (i, c) in enumerate(clusterids(exp))
        exclude[i] = !(c in idset)
    end
    deleteat!(exp.clusters, exclude)
    deleteat!(exp.clusterids, exclude)
end

