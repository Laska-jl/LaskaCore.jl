##############################################################
#
# Functions for converting times (ms <-> samprate for example)
#
##############################################################

# sample frequency/time conversion

const TUnit = Quantity{<:Number,Unitful.𝐓}

"""
    timetosamplerate(cluster::T, time::U) where {T<:AbstractCluster, U<:Quantity{<:Number, Unitful.𝐓}}

Convert a `time` of some unit from `Unitful.jl` into the samplerate of the `cluster`. For example, converting 10ms to a samplerate of 30 000Hz will yield a value of 300.0.
The returned range will not be of any `Unitful` unit.

# Examples

```julia
using LaskaCore
using Unitful

c = getcluster(exp, 33) # Get cluster '33' from an AbstractExperiment

timetosamplerate(c, 10u"ms") # Convert 10ms into the samplerate of the cluster
```
"""
function timetosamplerate(
    cluster::T,
    time::U,
) where {T<:AbstractCluster,U<:TUnit}
    samp::Float64 = @views info(cluster, "samprate")
    samp * ustrip(u"s", time)
end

"""
    timetosamplerate(cluster::T, time::U) where {T<:AbstractCluster,U<:AbstractRange{<:Quantity{<:Number,Unitful.𝐓}}}

Convert a `time`-range of some unit from `Unitful.jl` into the samplerate of the `cluster`. For example, converting 10ms to a samplerate of 30 000Hz will yield a value of 300.0.
The returned range will not be of any `Unitful` unit.

# Examples

```julia
using LaskaCore
using Unitful

c = getcluster(exp, 33) # Get cluster '33' from an AbstractExperiment

timetosamplerate(c, (0:10)u"ms") # Convert 0:10ms into the samplerate of the cluster
```
"""
function timetosamplerate(
    cluster::T,
    time::U,
) where {T<:AbstractCluster,U<:AbstractRange{<:TUnit{<:Number}}}
    samp::Float64 = @views info(cluster, "samprate")
    (samp*ustrip(u"s", time[begin])):(samp*ustrip(u"s", time[end]))
end

function timetosamplerate(
    V::T,
    time::U,
) where {T<:AbstractSpikeVector,U<:AbstractRange{<:TUnit{<:Number}}}
    samp = @views samplerate(V)
    return samp*ustrip(u"s", time[begin]):samp*ustrip(u"s", time[end])
end

function timetosamplerate(
    V::T,
    time::U,
) where {T<:AbstractSpikeVector,U<:TUnit{<:Number}}
    samp = @views samplerate(V)
    return samp * ustrip(u"s", time[begin])
end

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

In-place function for converting a vector to a samplerate. Multiplies each `in` with `factor` and puts the result in `out`.

Only for internal use.
"""
function mstosamplerate!(
    out::Vector{Float64},
    in::Vector{T},
    factor::Float64,
) where {T<:Real}
    @assert length(out) == length(in)
    @inbounds for i in eachindex(out)
        out[i] = in[i] * factor
    end
    return out
end




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

function sampleratetoms!(out::Vector{T}, factor::T) where {T<:AbstractFloat}
    @inbounds for n in eachindex(out)
        out[n] *= factor
    end
end
