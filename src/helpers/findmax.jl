

"""
    maxval(vec::Vector{Vector{T}}, init::T=0) where {T<:Real}

Find the maximum value in a Vector of Vectors.

"""
function maxval(vec::Vector{Vector{T}}, init::T=0) where {T<:Real}
    max::T = init
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
    minval(vec::Vector{Vector{T}}, init::T=0) where {T<:Real}

Find the minimum value in a Vector of Vectors.

"""
function minval(vec::Vector{Vector{T}}, init::T=0) where {T<:Real}
    min::T = init
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
