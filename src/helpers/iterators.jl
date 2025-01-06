#=
Iteration utilities for spikevectors, clusters etc.
=#

"""
    SpikeTrainIterator(
        spikes::SpikeVector{T},
        t_start::T,
        t_end::T,
    ) where {T<:Integer}

Iterator which lazily iterates over all sequential times between `t_start` and `t_end`, returning 1 at times where there is a spike.
Currently only supports `Integer` time steps. A `SpikeVector` with spikes at t = 2, 3, 5 with `t_start` = 0 and `t_end` = 5 would yield:            
`[0, 0, 1, 1, 0, 1]`

# Example

```julia
it = SpikeTrainIterator(spikes, 0, 1000)

for s in it
    # Do something based on if s = 1 or s = 0
end

```
"""
struct SpikeTrainIterator{T<:Integer}
    spikes::SpikeVector{T}
    t_start::T
    t_end::T

    function SpikeTrainIterator(
        spikes::SpikeVector{T},
        t_start::T,
        t_end::T,
    ) where {T<:Integer}
        if t_start > spikes[begin] || t_end < spikes[end]
            return new{T}(
                LaskaCore.spikes_in_timerange(spikes, t_start, t_end),
                t_start,
                t_end,
            )
        else
            return new{T}(spikes, t_start, t_end)
        end
    end
end

function Base.iterate(v::SpikeTrainIterator)
    # If start time is the same as the first spike, return 1 and a tuple (1, (state, second spike)) where 1 = there was a spike there
    v.t_start == v.spikes[begin] && return (1, (2, 2))
    # ...else return (0, (state, first spike)) where 0 = there was not a spike there
    return (0, (2, 1))
end

function Base.iterate(v::SpikeTrainIterator, state)
    state[1] > v.t_end && return nothing
    state[1] > v.spikes[end] && return (0, (state[1] + 1, state[2]))
    state[1] == v.spikes[state[2]] && return (1, (state[1] + 1, state[2] + 1))
    return (0, (state[1] + 1, state[2]))
end

function Base.eltype(::SpikeTrainIterator{T}) where {T}
    return T
end

function Base.length(v::SpikeTrainIterator)
    v.t_end - v.t_start
end

function Base.size(v::SpikeTrainIterator)
    (v.t_end - v.t_start,)
end
