###############################################
#
# Get a vector of all spikes at a certain depth
#
###############################################

"""
    spikesatdepth(p::PhyOutput{T}, depth) where {T<:Real}
    spikesatdepth(p::PhyOutput{T}, depth::NTuple{2}) where {T<:Real}
    spikesatdepth(p::PhyOutput{T}, depth::Set) where {T<:Real}


Returns a `Vector{T}` of all spiketimes at/in `depth`. Spiketimes are **not sorted** on return.

The included depths are controlled by the type of the `depth` variable:

- A **single number** returns only the spikes of clusters at that exact depth.
- A **Tuple** with 2 entries returns all clusters at depths between (and including) the values.
- A **Set** returns the spiketimes of clusters with the exact depths in the Set.

Returns an empty `SpikeVector` if no clusters match the specified depth(s).
"""
function spikesatdepth(p::PhyOutput{T}, depth) where {T<:Real}
    out::Vector{T} = T[]
    for cluster in clustervector(p)
        if info(cluster, "depth") == depth
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end

function spikesatdepth(p::PhyOutput{T}, depth::NTuple{2}) where {T<:Real}
    out::Vector{T} = T[]
    for cluster in clustervector(p)
        if depth[1] <= info(cluster, "depth") <= depth[2]
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end


function spikesatdepth(p::PhyOutput{T}, depth::Set) where {T<:Real}
    out::Vector{T} = T[]
    for cluster in clustervector(p)
        if info(cluster, "depth") in depth
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end

# Versions for relativespikes

"""
    spikesatdepth(p::RelativeSpikes{T}, depth) where {T<:Real}
    spikesatdepth(p::RelativeSpikes{T}, depth::NTuple{2}) where {T<:Real}
    spikesatdepth(p::RelativeSpikes{T}, depth::Set) where {T<:Real}


Returns a single `RelativeSpikes{T,U}`-vector containing all spiketimes at/in `depth`. Spiketimes are **not sorted** on return.
Samplerates of all clusters in `p` must be equal.

The included depths are controlled by the type of the `depth` variable:

- A **single number** returns only the spikes of clusters at that exact depth.
- A **Tuple** with 2 entries returns all clusters at depths between (and including) the values.
- A **Set** returns the spiketimes of clusters with the exact depths in the Set.

Returns an empty `RelativeSpikeVector` if no clusters match the specified depth(s).
"""
function spikesatdepth(p::RelativeSpikes{T}, depth) where {T<:Real}
    samplerates = LaskaCore.samplerates(p)
    samp = all(samplerates .== samplerates[1]) ? samplerates[1] : throw(ErrorException("Samplerates of all clusters not equal"))
    out = RelativeSpikeVector([T[] for _ in 1:ntrigs(p)], samp)
    for cluster in clustervector(p)
        if info(cluster, "depth") == depth
            times = spiketimes(cluster)
            for i in eachindex(spiketimes(cluster))
                push!(out[i], times[i]...)
            end
        end
    end
    return out
end

function spikesatdepth(p::RelativeSpikes{T}, depth::NTuple{2}) where {T<:Real}
    samplerates = LaskaCore.samplerates(p)
    samp = all(samplerates .== samplerates[1]) ? samplerates[1] : throw(ErrorException("Samplerates of all clusters not equal"))
    out = RelativeSpikeVector([T[] for _ in 1:ntrigs(p)], samp)
    for cluster in clustervector(p)
        if depth[1] <= info(cluster, "depth") < depth[2]
            times = spiketimes(cluster)
            for i in eachindex(spiketimes(cluster))
                if length(times[i]) > 0
                    push!(out[i], times[i]...)
                end
            end
        end
    end
    return out
end


function spikesatdepth(p::RelativeSpikes{T}, depth::Set) where {T<:Real}
    samplerates = LaskaCore.samplerates(p)
    samp = all(samplerates .== samplerates[1]) ? samplerates[1] : throw(ErrorException("Samplerates of all clusters not equal"))
    out = RelativeSpikeVector([T[] for _ in 1:ntrigs(p)], samp)
    for cluster in clustervector(p)
        if info(cluster, "depth") in depth
            times = spiketimes(cluster)
            for i in eachindex(spiketimes(cluster))
                if length(times[i]) > 0
                    push!(out[i], times[i]...)
                end
            end
        end
    end
    return out
end
