##############################################################
#
# Functions for converting times (ms <-> samprate for example)
#
##############################################################


"""
    mstosamplerate(ms::Int64, experiment::T) where {T<:AbstractExperiment}
    mstosamplerate(vec::Vector{T}, samplerate::U) where {T<:Real,U<:Real}
    mstosamplerate(ms::Int64, samplerate::T) where {T<:Real}

Convert a value in ms to samplerate of `experiment`. A custom samplerate may be provided instead of an `experiment`

"""
function mstosamplerate(ms::Int64, experiment::T) where {T<:AbstractExperiment}
    return ms * parse(Float64, experiment.meta["imSampRate"]) * 0.001
end

function mstosamplerate(ms::Int64, samplerate::T) where {T<:Real}
    return ms * samplerate * 0.001
end
# Version for vectors

function mstosamplerate(vec::Vector{T}, samplerate::U) where {T<:Real,U<:Real}
    fac = samplerate * 0.001
    out::Vector{Float64} = Vector{Float64}(undef, length(vec))
    mstosamplerate!(out, vec, fac)
    return out
end

# NOTE: Replace with standard scalar?
"""
    mstosamplerate!(out::Vector{Float64}, in::Vector{T}, factor::Float64) where {T<:Real}

In-place function for converting a vector to a samplerate. Very simply multiplies each `in` with `factor` and puts the result in `out`.

Only for internal use.
"""
function mstosamplerate!(out::Vector{Float64}, in::Vector{T}, factor::Float64) where {T<:Real}
    @assert length(out) == length(in)
    @inbounds for i in eachindex(out)
        out[i] = in[i] * factor
    end
    return out
end


# FIX: reverse argument order


"""
    sampleratetoms(ms::Int64, experiment::T) where {T<:AbstractExperiment}
    sampleratetoms(ms::Int64, samplerate::T) where {T<:Real}
    sampleratetoms(vec::Vector{T}, samplerate::U) where {T<:Real,U<:Real}

Convert a time measured in samplerate of an `AbstractExperiment` to ms. A direct `samplerate` may be provided instead of an `AbstractExperiment`.
"""
function sampleratetoms(ms::Int64, experiment::T) where {T<:AbstractExperiment}
    return ms / (experiment.meta["imSampRate"] * 0.001)
end

function sampleratetoms(ms::Int64, samplerate::T) where {T<:Real}
    return ms / (samplerate * 0.001)
end

function sampleratetoms(vec::Vector{T}, samplerate::U) where {T<:Real,U<:Real}
    fac = 1 / (samplerate * 0.001)
    out::Vector{Float64} = deepcopy(vec)
    sampleratetoms!(out, fac)
    return out
end

function sampleratetoms!(out::Vector{Float64}, factor::Float64)
    @inbounds for n in eachindex(out)
        out[n] *= factor
    end
end

