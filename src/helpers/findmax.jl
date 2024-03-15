"""
    maxval(vec::RelativeSpikeVector, init::T=0) where {T<:Real}

Find the maximum value in a [`RelativeSpikeVector`](@ref).

"""
function maxval(vec::RelativeSpikeVector, init::T=0) where {T}
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

function maxval(vec::Vector{Vector{T}}, init::U=0) where {T,U}
    max::U = init
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
function minval(vec::RelativeSpikeVector, init::T=0) where {T}
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

function minval(vec::Vector{Vector{T}}, init::U=0) where {T,U}
    min = init
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
