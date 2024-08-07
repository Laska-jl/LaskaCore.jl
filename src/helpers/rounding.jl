###############
#
# Round up/down
#
###############

"""
    roundup(value::T, interval::N) where T<:Real where N<:Real

Rounds `value` up to the nearest greater multiple of `interval`.

"""
function roundup(value::T, interval::N) where {T,N}
    return N(ceil(value / interval) * interval)
end

"""

    rounddown(value::T, interval::N) where {T<:Real} where {N<:Real}

Rounds `value` down to the nearest lesser multiple of `interval`.

"""
function rounddown(value::T, interval::N) where {T,N}
    return N(floor(value / interval) * interval)
end

"""

    arbitraryround(value::T, interval::N) where {T<:Real} where {N<:Real}

Rounds `value` to the nearest multiple of `interval`.

"""
function arbitraryround(value::T, interval::N) where {T,N}
    return N(Base.round(value / interval) * interval)
end


function roundupvec(in::Vector{T}, period::N) where {T,N}
    out = Vector{N}(undef, length(in))
    roundupvec!(out, in, period)
    return out
end

function roundupvec!(out::Vector{T}, in::Vector{N}, period::T) where {T,N}
    for n in eachindex(out)
        out[n] = roundup(in[n], period)
    end
end

function rounddownvec(in::Vector{T}, period::N) where {T,N}
    out = Vector{N}(undef, length(in))
    rounddownvec!(out, in, period)
    return out
end

function rounddownvec!(out::Vector{T}, in::Vector{N}, period::T) where {T,N}
    for n in eachindex(out)
        out[n] = rounddown(in[n], period)
    end
end

function arbitraryroundvec(in::Vector{T}, period::N) where {T,N}
    out = Vector{N}(undef, length(in))
    arbitraryroundvec!(out, in, period)
    return out
end

function arbitraryroundvec!(out::Vector{T}, in::Vector{N}, period::T) where {T,N}
    @inbounds for n in eachindex(out)
        out[n] = arbitraryround(in[n], period)
    end
end


