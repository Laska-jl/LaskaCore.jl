"""
    maxval(vec::RelativeSpikeVector, init::T=0) where {T<:Real}

Find the maximum value in a [`RelativeSpikeVector`](@ref).

"""
function maxval(vec::RelativeSpikeVector)
    max = vec[1][1]
    @inbounds @views begin
        for n in eachindex(vec)
            for i in eachindex(vec[n])
                if vec[n][i] > max
                    max = vec[n][i]
                end
            end
        end
    end
    return max
end

function maxval(vec::Vector{Vector{T}}) where {T}
    max = vec[1][1]
    @inbounds @views begin
        for n in eachindex(vec)
            for i in eachindex(vec[n])
                if vec[n][i] > max
                    max = vec[n][i]
                end
            end
        end
    end
    return max
end

"""
    minval(vec::RelativeSpikeVector, init::T=0) where {T<:Real}

Find the minimum value in a [`RelativeSpikeVector`](@ref).

"""
function minval(vec::RelativeSpikeVector)
    min = vec[1][1]
    @inbounds @views begin
        for n in eachindex(vec)
            for i in eachindex(vec[n])
                if vec[n][i] < min
                    min = vec[n][i]
                end
            end
        end
    end
    return min
end

function minval(vec::Vector{Vector{T}}) where {T}
    min = vec[1][1]
    @inbounds @views begin
        for n in eachindex(vec)
            for i in eachindex(vec[n])
                if vec[n][i] < min
                    min = vec[n][i]
                end
            end
        end
    end
    return min
end

"""
    function extremevals(vec::Vector{Vector{T}}) where {T,U}
    function extremevals(vec::RelativeSpikeVector)

Returns the extrema of a [`LaskaCore.RelativeSpikeVector`](@ref) or a `Vector{Vector{T}}`, similar to the `extrema` function from `Base`.
"""
function extremevals end

function extremevals(vec::Vector{Vector{T}}) where {T}
    return minval(vec), maxval(vec)
end

function extremevals(vec::RelativeSpikeVector)
    return minval(vec), maxval(vec)
end
