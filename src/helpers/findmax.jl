"""
    maxval(vec::RelativeSpikeVector, init::T=0) where {T<:Real}

Find the maximum value in a [`RelativeSpikeVector`](@ref).

"""
function maxval(vec::RelativeSpikeVector{T}) where {T}
    maximum(maximum.(vec))
end

function maxval(vec::Vector{Vector{T}}) where {T}
    maximum(maximum.(vec))
end

"""
    minval(vec::RelativeSpikeVector, init::T=0) where {T<:Real}

Find the minimum value in a [`RelativeSpikeVector`](@ref).

"""
function minval(vec::RelativeSpikeVector{T}) where {T}
    minimum(minimum.(vec))
end

function minval(vec::Vector{Vector{T}}, init::T=zero(T)) where {T}
    minimum(minimum.(vec))
end

"""
    function extremevals(vec::Vector{Vector{T}}, init::U=0) where {T,U}
    function extremevals(vec::RelativeSpikeVector, init::U=0) where {U}

Returns the extrema of a [`LaskaCore.RelativeSpikeVector`](@ref) or a `Vector{Vector{T}}`, similar to the `extrema` function from `Base`.
"""
function extremevals end

function extremevals(vec::Vector{Vector{T}}) where {T}
    return minval(vec), maxval(vec)
end

function extremevals(vec::RelativeSpikeVector{T}) where {T}
    return minval(vec), maxval(vec)
end
