###############
#
# Round up/down
#
###############

"""

    roundup(value::T, interval::N) where T<:Real where N<:Real

Rounds `value` up to the nearest greater multiple of `interval`.


"""
function roundup(value::T, interval::N) where {T<:Real} where {N<:Real}
    return N(ceil(value / interval) * interval)
end

"""

    rounddown(value::T, interval::N) where {T<:Real} where {N<:Real}

Rounds `value` down to the nearest lesser multiple of `interval`.

"""
function rounddown(value::T, interval::N) where {T<:Real} where {N<:Real}
    return N(floor(value / interval) * interval)
end

"""

    arbitraryround(value::T, interval::N) where {T<:Real} where {N<:Real}

Rounds `value` to the nearest multiple of `interval`.

"""
function arbitraryround(value::T, interval::N) where {T<:Real} where {N<:Real}
    return N(Base.round(value / interval) * interval)
end


function roundupvec(in::Vector{N}, period::T) where {T<:Real} where {N<:Real}
    out = Vector{T}(undef, length(in))
    roundupvec!(out, in, period)
    return out
end

function roundupvec!(out::Vector{T}, in::Vector{N}, period::T) where {N<:Real} where {T<:Real}
    @inbounds for n in eachindex(out)
        out[n] = roundup(in[n], period)
    end
end

function rounddownvec(in::Vector{N}, period::T) where {T<:Real} where {N<:Real}
    out = Vector{T}(undef, length(in))
    rounddownvec!(out, in, period)
    return out
end

function rounddownvec!(out::Vector{T}, in::Vector{N}, period::T) where {N<:Real} where {T<:Real}
    @inbounds for n in eachindex(out)
        out[n] = rounddown(in[n], period)
    end
end

function arbitraryroundvec(in::Vector{N}, period::T) where {T<:Real} where {N<:Real}
    out = Vector{T}(undef, length(in))
    arbitraryroundvec!(out, in, period)
    return out
end

function arbitraryroundvec!(out::Vector{T}, in::Vector{N}, period::T) where {N<:Real} where {T<:Real}
    @inbounds for n in eachindex(out)
        out[n] = arbitraryround(in[n], period)
    end
end


