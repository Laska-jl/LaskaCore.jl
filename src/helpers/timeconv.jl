##############################################################
#
using CSV: Error
# Functions for converting times (ms <-> samprate for example)
#
##############################################################

# sample frequency/time conversion

# Convenience type for Unitful units of time
const TUnit{T} = Quantity{T,Unitful.𝐓}
const FreeTUnit{T} = Unitful.FreeUnits{T,Unitful.𝐓}

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
) where {T<:Union{<:AbstractCluster,<:AbstractSpikeVector},U<:TUnit}
    samplerate(cluster) * ustrip(Float64, u"s", time)
end

"""
    timetosamplerate(cluster::T, time::U) where {T<:AbstractCluster,U<:AbstractRange{<:Quantity{<:Number,Unitful.𝐓}}}

Convert a `time`-range of some unit from `Unitful.jl` into the samplerate of the `cluster`. For example, converting 10ms to a samplerate of 30 000Hz will yield a value of 300.0.
The returned "time" will not be of any `Unitful` unit.

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
) where {T<:Union{<:AbstractCluster,<:AbstractSpikeVector},U<:AbstractRange{<:TUnit{<:Number}}}
    samp = samplerate(cluster)
    (samp*ustrip(Float64, u"s", time[begin])):(samp*ustrip(Float64, u"s", time[end]))
end


function timetosamplerate(V::T, time::U) where {T<:Union{<:AbstractCluster,<:AbstractSpikeVector},U<:StepRange{<:TUnit{<:Number}}}
    samp = samplerate(V)
    lower = samp * ustrip(Float64, u"s", time[begin])
    step = samp * ustrip(Float64, u"s", time.step)
    upper = samp * ustrip(Float64, u"s", time[end])
    return lower:step:upper
end

# Version for non-Int stepranges with Unitful units.
function timetosamplerate(V::T, time::U) where {T<:Union{<:AbstractCluster,<:AbstractSpikeVector},U<:StepRangeLen{<:TUnit{<:Number}}}
    samp = samplerate(V)
    lower = samp * ustrip(Float64, u"s", time[begin])
    step = samp * ustrip(Float64, u"s", time.step)
    upper = samp * ustrip(Float64, u"s", time[end])
    return lower:step:upper
end

# Version for entire experiment, returns times for each cluster
"""
    timetosampleratevec(p::AbstractExperiment, time)

Converts `time` to the sample rate of each cluster in `p` and returns the result as a Vector.

"""
function timetosampleratevec(p::AbstractExperiment, time)
    clusters = clustervector(p)
    [timetosamplerate(clusters[i], time) for i in eachindex(clusters)]
end

# Version which returns a single time for an experiment if all times are equal
"""
    timetoequalsamplerate(p::AbstractExperiment, time::TUnit)

Converts `time` into the samplerate of `p` if all clusters in `p` have the same samplerate.
Otherwise returns an error.
"""
function timetosamplerate(p::AbstractExperiment, time)
    all = timetosampleratevec(p, time)
    if length(unique(all)) != 1
        throw(ErrorException("Sampletimes of all clusters not equal"))
    end
    return all[1]
end

"""
    sampleratetotime(samplerate::T, time::T, unit::U) where {T<:Real, U<:LaskaCore.FreeTUnit}

Convert a time represented in a samplerate to another unit. `unit` should ba `Unitful.jl` unit of time. Always returns a normal `Float64`.
`time' may be a single number or a `Vector` of such.


# Example

```julia
# Convert time = 450000 sampled at 30000 Hz to ms (i.e 1500)
LaskaCore.sampleratetotime(30000, 45000, u"ms")

```
"""
function sampleratetotime(samplerate::Real, time::Real, unit::U) where {U<:LaskaCore.FreeTUnit}
    ustrip(Float64, unit, (time / samplerate)u"s")
end

function sampleratetotime!(out::Vector{<:Real}, samplerate::Real, times::T, unit::U) where {T<:AbstractVector{<:Real},U<:LaskaCore.FreeTUnit}
    for i in eachindex(out)
        out[i] = ustrip(Float64, unit, (time[i] / samplerate)u"s")
    end
end

function sampleratetotime(samplerate::Real, time::T, unit::U) where {T<:AbstractVector{<:Real},U<:LaskaCore.FreeTUnit}
    out = similar(time, Float64)
    @simd for i in eachindex(out)
        @views out[i] = sampleratetotime(samplerate, time[i], unit)
    end
    return out
end

