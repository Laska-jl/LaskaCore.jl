#=
Functions for filtering SpikeVectors
=#

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

function filterduplicatespikes!(v::SpikeVector)
    inds = zeros(UInt64, 0)
    v[1] == v[2] && push!(inds, 1)
    for i = 2:length(v)-1
        v[i] == v[i+1] && push!(inds, i)
    end
    v[end-1] == v[end] && push!(inds, length(v))
    deleteat!(v, inds)
end
