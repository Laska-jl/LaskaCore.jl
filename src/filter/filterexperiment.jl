# Functions for filtering experiments

"""
    Base.filter!((var, f)::Pair, exp::AbstractExperiment)

Filter the `Cluster`s of `exp` based off a column in its info `DataFrame`. The syntax is the same as when filtering a "normal" DataFrame.

# Example

In order to remove all clusters of `exp` with a depth less than 1000:

```julia
filter!(:depth => x -> x >= 1000.0, exp)
```
"""
function Base.filter!((var, f)::Pair, exp::AbstractExperiment)
    allids = deepcopy(clusterids(exp))
    filt = [true for _ in eachindex(allids)]
    DataFrames.filter!(var => f, info(exp))
    valid = Set(info(exp, :cluster_id))
    for i in eachindex(allids)
        if !(allids[i] in valid)
            filt[i] = false
        end
    end
    exp.clusters = clustervector(exp)[filt]
end

