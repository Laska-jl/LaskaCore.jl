#=
Functions for filtering SpikeVectors
=#

"""
    filterduplicatespikes(v::SpikeVector)

Finds all spikes with identical spiketimes in `v` and removes one of them. Returns a new `SpikeVector`, leaving `v` unchanged.
For mutating version see [`LaskaCore.filterduplicatespikes!`](@ref).

"""
function filterduplicatespikes(v::SpikeVector)
    inds = zeros(UInt64, 0)
    v[1] == v[2] && push!(inds, 1)
    for i = 2:length(v)-1
        v[i] == v[i+1] && push!(inds, i)
    end
    v[end-1] == v[end] && push!(inds, length(v))
    out = deepcopy(v)
    deleteat!(out, inds)
    return out
end

"""
    filterduplicatespikes!(v::SpikeVector)

Finds all spikes with identical spiketimes in `v` and removes one of them. Mutates `v`, for non-mutating version see [`LaskaCore.filterduplicatespikes`](@ref).

"""
function filterduplicatespikes!(v::SpikeVector)
    inds = zeros(UInt64, 0)
    v[1] == v[2] && push!(inds, 1)
    for i = 2:length(v)-1
        v[i] == v[i+1] && push!(inds, i)
    end
    v[end-1] == v[end] && push!(inds, length(v))
    deleteat!(v, inds)
end



"""
    filterstims(spikes::AbstractVector{T}, stimtimes::AbstractVector{T}, window::T)

- `spikes`: Vector of spiketimes
- `stimtimes`: Vector of stimulation times
- `window`: Time before and after each stimtime to remove
"""
function filterstims(spikes::AbstractVector{T}, stimtimes::AbstractVector{T}, window::T) where {T}
	stimwindows = [stimtimes[i] - window:stimtimes[i]+window for i in eachindex(stimtimes)]

    out = Vector{T}[]
    startind = 1
    for i in eachindex(stimwindows)
        cur_view = @view spikes[startind:end]
        inds = findall(x -> x < stimwindows[i][end], cur_view)
        v = @view cur_view[inds]
        push!(out, filter(x -> x < stimwindows[i][begin], v))
        startind += inds[end]
    end
    return out
end

function anyin(v)
    x -> any(in.(x, v))
end
