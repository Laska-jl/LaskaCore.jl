#----------------------------
# Iterate over spiketimes etc
#----------------------------

# Iterate over spiketimes of Cluster

function Base.iterate(cluster::T) where {T<:AbstractCluster}
    spikes = spiketimes(cluster)
    state = [1, 1, length(spikes)]
    if state[1] < spikes[state[2]]
        state[1] += 1
        return 0, state
    else
        state[1] += 1
        state[2] += 1
        return state[3] < state[1] ? (1, state) : nothing
    end
end

function Base.iterate(cluster::T, state::Vector{Int}) where {T<:AbstractCluster}
    spikes = spiketimes(cluster)
    if state[1] < spikes[state[2]]
        state[1] += 1
        return 0, state
    else
        state[1] += 1
        state[2] += 1
        return state[3] < state[1] ? (1, state) : nothing

    end

end
