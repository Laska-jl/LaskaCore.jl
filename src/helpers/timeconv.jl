##############################################################
#
# Functions for converting times (ms <-> samprate for example)
#
##############################################################

# sample frequency/time conversion

# Convenience type for Unitful units of time
const TUnit{T} = Quantity{T,Unitful.𝐓}

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


