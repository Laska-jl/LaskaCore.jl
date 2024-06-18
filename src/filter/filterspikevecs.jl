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
