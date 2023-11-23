###############################################
#
# Get a vector of all spikes at a certain depth
#
###############################################

"""

    spikesatdepth(experiment::E{T}, depth::N) where {E<:AbstractExperiment,T<:Real,N<:Real}
    spikesatdepth(experiment::E{T}, depth::Tuple{2,N}) where {E<:AbstractExperiment,T<:Real,N<:Real}
    spikesatdepth(experiment::E{T}, depth::Set{N}) where {E<:AbstractExperiment,T<:Real,N<:Real}

Returns a `Vector{T}` of all spiketimes at/in `depth`.

The included depths are controlled by the type of the `depth` variable:                 

- A **single number** returns only the spikes of clusters at that exact depth.                  
- A **Tuple** with 2 entries returns all clusters at depths between (and including) the values.                  
- A **Set** returns the clusters with the exact depths in the Set.
"""
function spikesatdepth(p::PhyOutput{T}, depth::N) where {T<:Real} where {N<:Real}
    out::Vector{T} = T[]
    for cluster in clustervector(p)
        if parse(N, info(cluster, "depth")) == depth
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end

function spikesatdepth(p::PhyOutput{T}, depth::NTuple{2,N}) where {T<:Real} where {N<:Real}
    out::Vector{T} = T[]
    for cluster in clustervector(p)
        if depth[1] <= parse(N, info(cluster, "depth")) <= depth[2]
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end


function spikesatdepth(p::PhyOutput{T}, depth::Set{N}) where {T<:Real} where {N<:Real}
    out::Vector{T} = T[]
    for cluster in clustervector(p)
        if parse(N, info(cluster, "depth")) in depth
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end

# Versions for relativespikes

function spikesatdepth(p::RelativeSpikes{T}, depth::N) where {T<:Real} where {N<:Real}
    out::Vector{Vector{T}} = Vector{Vector{T}}(undef, 0)
    for cluster in clustervector(p)
        if parse(N, info(cluster, "depth")) == depth
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end

function spikesatdepth(p::RelativeSpikes{T}, depth::NTuple{2,N}) where {T<:Real} where {N<:Real}
    out::Vector{Vector{T}} = Vector{Vector{T}}(undef, 0)
    for cluster in clustervector(p)
        if depth[1] <= parse(N, info(cluster, "depth")) < depth[2]
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end


function spikesatdepth(p::RelativeSpikes{T}, depth::Set{N}) where {T<:Real} where {N<:Real}
    out::Vector{Vector{T}} = Vector{Vector{T}}(undef, 0)
    for cluster in clustervector(p)
        if parse(N, info(cluster, "depth")) in depth
            out = vcat(out, spiketimes(cluster))
        end
    end
    return out
end
