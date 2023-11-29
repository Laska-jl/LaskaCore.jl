


"""

    normalize(vec::Vector{T}) where {T<:Real}
    normalize(vec::Vector{T}, min::U, max::U) where {T<:Real,U<:Real}
    normalize(vec::Vector{T}, min::U, max::U, minmax::NTuple{2,U}) where {T<:Real,U<:Real}
    normalize(vec::Vector{T}, minmax::NTuple{2,U}) where {T<:Real,U<:Real}


Returns a min-max normalized version of `vec`. Defaults to 0-1 which can be customized with (`min`/`max`).
Custom values replacing the actual min/max values in the vector may also be supplied in a Tuple with exactly 2 entries.

See [`LaskaCore.normalize!`](@ref) for an in-place version
"""
function normalize(vec::Vector{T}) where {T<:Real}
    out::Vector{Float64} = deepcopy(vec)
    normalize!(out)
    return out
end

function normalize(vec::Vector{T}, min::U, max::U) where {T<:Real,U<:Real}
    out::Vector{Float64} = deepcopy(vec)
    normalize!(out, min, max)
    return out
end

function normalize(vec::Vector{T}, min::U, max::U, minmax::NTuple{2,U}) where {T<:Real,U<:Real}
    out::Vector{Float64} = deepcopy(vec)
    normalize!(out, min, max, minmax)
    return out
end

function normalize(vec::Vector{T}, minmax::NTuple{2,U}) where {T<:Real,U<:Real}
    out::Vector{Float64} = deepcopy(vec)
    normalize!(out, minmax)
    return out
end

"""

Normalize an `AbstractFloat` vector in place.

"""
function normalize!(vec::Vector{T}) where {T<:AbstractFloat}
    minx = minimum(vec)
    diffx = maximum(vec) - minx
    @inbounds for v in eachindex(vec)
        vec[v] = ((vec[v] - minx)) / diffx
    end
end

function normalize!(vec::Vector{T}, min::U, max::U) where {T<:AbstractFloat,U<:Real}
    minx = minimum(vec)
    diffx = maximum(vec) - minx
    rangediff = max - min
    @inbounds for v in eachindex(vec)
        vec[v] = min + ((vec[v] - minx) * rangediff) / diffx
    end
end

function normalize!(vec::Vector{T}, min::U, max::U, minmax::NTuple{2,U}) where {T<:AbstractFloat,U<:Real}
    minx = minmax[1]
    diffx = minmax[2] - minmax[1]
    rangediff = max - min
    @inbounds for v in eachindex(vec)
        vec[v] = min + ((vec[v] - minx) * rangediff) / diffx
    end
end

function normalize!(vec::Vector{T}, minmax::NTuple{2,U}) where {T<:AbstractFloat,U<:Real}
    diffx = minmax[2] - minmax[1]
    @inbounds for v in eachindex(vec)
        vec[v] = ((vec[v] - minmax[1])) / diffx
    end
end




