#=========================================================================
Vector type for holding spiketimes and the Hz/unit of their representation
=========================================================================#

const HzType = Quantity{<:Number,Unitful.𝐓^-1}

"""
    abstract type AbstractSpikeVector{T,U} <: AbstractArray{T,1} end

Supertype for `Vector`s holding spiketimes and the samplerate at which they were collected/are represented as.
"""
abstract type AbstractSpikeVector{T,U} <: AbstractArray{T,1} end

"""
    struct SpikeVector{T,U} <: AbstractSpikeVector{T,U}
        spiketimes::Vector{T}
        samplerate::U
    end

Struct for holding spiketimes and the sample rate in Hz.

Since this is not a mutable struct the sample rate cannot be changed in-place. The spiketimes themselves may however be changed since they are stored in a mutable `Vector{T}`.
For conversion to a new sample rate, creating a new `RelativeSpikeVector` using [`LaskaCore.convertsamplerate`](@ref) is recommended.
"""
struct SpikeVector{T,U} <: AbstractSpikeVector{T,U}
    spiketimes::Vector{T}
    samplerate::U
end
export SpikeVector

"""
    struct RelativeSpikeVector{T,U} <: AbstractSpikeVector{T,U}
        spiketimes::Vector{Vector{T}}
        samplerate::U
    end

Struct for holding spiketimes relative to some trigger as well as their samplerate.

The spiketimes themselves are held in a `Vector{Vector{T}}` (A vector of vectors) where the inner vector `n` holds spiketimes relative to the `n`:th trigger event.

Since this is not a mutable struct the sample rate cannot be changed in-place. The spiketimes themselves may however be changed since they are stored in a mutable `Vector{T}`.
For conversion to a new sample rate, creating a new instance of `RelativeSpikeVector` using [`LaskaCore.convertsamplerate`](@ref) is recommended.

# Indexing

Since this is essentially just a vector of vectors it can be indexed as `vec[i][j]` where `i` refers to the `i`:th trigger events and `j` to the `j`:th spiketime tied to that trigger event.
"""
struct RelativeSpikeVector{T,U} <: AbstractSpikeVector{T,U}
    spiketimes::Vector{Vector{T}}
    samplerate::U
end
export RelativeSpikeVector

#===============================
Functions on AbstractSpikeVector
===============================#


"""
    Base.vec(V::AbstractSpikeVector{T,U}) where {T,U}

Returns a **view** into the spiketimes of `V` as a `Vector{T}`. Modifying the returned vector will therefore also modify `V`!

In order to modify the returned `Vector` without affecting `V` explicitly make a copy using `deepcopy(vec(V))`.
"""
@inline function Base.vec(V::AbstractSpikeVector{T,U}) where {T,U}
    Base.vec(V.spiketimes)
end

"""
    samplerate(V::AbstractSpikeVector)

Returns the sample rate of `V`.
"""
samplerate(V::AbstractSpikeVector) = V.samplerate

"""
    spiketimes(V::AbstractSpikeVector)

Returns the spiketimes of `V` as a `Vector`.
"""
spiketimes(V::AbstractSpikeVector) = V.spiketimes

@inline Base.size(V::AbstractSpikeVector) = Base.size(V.spiketimes)

@inline Base.getindex(V::AbstractSpikeVector, i::Int) = V.spiketimes[i]

Base.IndexStyle(::Type{<:AbstractSpikeVector}) = IndexLinear()

@inline Base.length(V::SpikeVector) = length(V.spiketimes)

@inline Base.append!(V::AbstractSpikeVector, items::AbstractVector) =
    Base.append!(V.spiketimes, items)

@inline Base.deleteat!(V::AbstractSpikeVector, i::Integer) = Base.deleteat!(V.spiketimes, i)
@inline Base.deleteat!(V::AbstractSpikeVector, r::AbstractUnitRange{<:Integer}) =
    Base.deleteat!(V.spiketimes, r)
@inline Base.deleteat!(V::AbstractSpikeVector, inds::AbstractVector{Bool}) = Base.deleteat!(V.spiketimes, inds)
@inline Base.deleteat!(V::AbstractSpikeVector, inds::AbstractVector) =
    Base.deleteat!(V.spiketimes, inds)
@inline Base.deleteat!(V::AbstractSpikeVector, inds) = Base.deleteat!(V.spiketimes, inds)

Base.firstindex(V::AbstractSpikeVector) = 1
Base.lastindex(V::AbstractSpikeVector) = length(V)

@inline Base.empty!(V::AbstractSpikeVector) = Base.empty!(V.spiketimes)

@inline function Base.insert!(V::AbstractSpikeVector{T}, i::Integer, item) where {T}
    Base.insert!(V.spiketimes, i, item)
end

@inline Base.keepat!(V::AbstractSpikeVector, m::AbstractVector{Bool}) =
    Base.keepat!(V.spiketimes, m)
@inline Base.keepat!(V::AbstractSpikeVector, inds) = Base.keepat!(V.spiketimes, inds)
@inline Base.keepat!(V::AbstractSpikeVector) = Base.pop!(V.spiketimes)

@inline Base.popat!(V::AbstractSpikeVector, i::Integer) = Base.popat!(V.spiketimes, i)
@inline Base.popat!(V::AbstractSpikeVector, i::Integer, default) =
    Base.popat!(V.spiketimes, i, default)

@inline Base.popfirst!(V::AbstractSpikeVector) = Base.popfirst!(V.spiketimes)

@inline Base.prepend!(V::AbstractSpikeVector, items::AbstractVector) =
    Base.prepend!(V.spiketimes, items)
@inline Base.prepend!(V::AbstractSpikeVector, iter) = Base.prepend!(V.spiketimes, iter)

@inline function Base.push!(V::AbstractSpikeVector{T}, item) where {T}
    Base.push!(V.spiketimes, item)
end

@inline function Base.pushfirst!(V::AbstractSpikeVector{T}, item) where {T}
    Base.push!(V.spiketimes, item)
end
@inline Base.pushfirst!(V::AbstractSpikeVector, iter...) =
    Base.pushfirst!(V.spiketimes, iter)

@inline Base.resize!(V::AbstractSpikeVector, nl::Integer) = Base.resize!(V.spiketimes, nl)

@inline Base.sizehint!(V::AbstractSpikeVector, sz::Integer) =
    Base.sizehint!(V.spiketimes, sz)

@inline Base.splice!(V::AbstractSpikeVector, i::Integer) = Base.splice!(V.spiketimes, i)
@inline Base.splice!(V::AbstractSpikeVector, i::Integer, ins) =
    Base.splice!(V.spiketimes, i, ins)
@inline Base.splice!(V::AbstractSpikeVector, r::AbstractUnitRange{<:Integer}) =
    Base.splice!(V.spiketimes, r)
@inline Base.splice!(V::AbstractSpikeVector, r::AbstractUnitRange{<:Integer}, ins) =
    Base.splice!(V.spiketimes, r, ins)
@inline Base.splice!(V::AbstractSpikeVector, inds) = Base.splice!(V.spiketimes, inds)

#=======================
Functions on SpikeVector
=======================#

"""
    (::Type{SpikeVector{T}})(::UndefInitializer, n, samplerate) where {T}


Initialize a `SpikeVector` of length `n` with `samplerate`
"""
function (::Type{SpikeVector{T}})(::UndefInitializer, n, samplerate) where {T}
    SpikeVector{T,typeof(samplerate)}(Vector{T}(undef, n), samplerate)
end

"""
    convertsamplerate(V::SpikeVector{T,U}, newsamplerate::N) where {T,U,N}

Convert `SpikeVector` `V` to a new samplerate. The type of the new vector will be `SpikeVector{Float64, typeof(newsamplerate)}`.

This function will create a new `SpikeVector` leaving `V` untouched.

# Example

In the example below, a `SpikeVector` with a sample rate of 30000Hz is converted to one with a sample rate of 1Hz.
```julia
julia> v = SpikeVector([i for i in 1:3000000], 30000)
3000000-element SpikeVector{Int64, Int64}:
       1
       2
       3
       4
       ⋮
 2999998
 2999999
 3000000

julia> LaskaCore.convertsamplerate(v, 1)
3000000-element SpikeVector{Float64, Int64}:
  30000.0
  60000.0
  90000.0
 120000.0
      ⋮
      8.999994e10
      8.999997e10
      9.0e10
```
"""
function convertsamplerate(V::SpikeVector{T,U}, newsamplerate) where {T,U}
    convfactor = samplerate(V) / newsamplerate
    outvec = Vector{Float64}(undef, length(V))
    @inbounds @views for i in eachindex(outvec)
        outvec[i] = V.spiketimes[i] * convfactor
    end
    return SpikeVector(outvec, newsamplerate)
end


"""
    convertsamplerate(V::SpikeVector{T,U}, newsamplerate, ::Type{N}) where {T,U,N}

If a specific output type of the spiketimes in the returned `SpikeVector` is desired, one may pass the type as the last argument.

This function will create a new `SpikeVector` leaving `V` untouched.

# Example

In the example below, a `SpikeVector{Float64, Int64}` with sample rate 30000Hz is converted to a `SpikeVector{Int64,Int64}` with sample rate 1Hz.
This requires that all resulting spiketimes of `V` after conversion is representable as `N`.

```julia
# Create a SpikeVector with spiketimes as Float64.
julia> v = SpikeVector([i for i in 1.0:30000.0], 30000)
30000-element SpikeVector{Float64, Int64}:
     1.0
     2.0
     3.0
     4.0
     ⋮
 29998.0
 29999.0
 30000.0

# Since the conversion becomes time * (30_000 / 1) all times will be representable as Int64:s
julia> LaskaCore.convertsamplerate(v, 1, Int64)
30000-element SpikeVector{Int64, Int64}:
     30000
     60000
     90000
    120000
         ⋮
 899940000
 899970000
 900000000
```
"""
function convertsamplerate(V::SpikeVector{T,U}, newsamplerate, ::Type{N}) where {T,U,N}
    convfactor = N(samplerate(V) / newsamplerate)
    outvec = Vector{N}(undef, length(V))
    @inbounds @views for i in eachindex(outvec)
        outvec[i] = V.spiketimes[i] * convfactor
    end
    return SpikeVector(outvec, newsamplerate)
end

# Interface implementation

@inline function Base.setindex!(V::SpikeVector{T}, val::T, i::Int) where {T}
    V.spiketimes[i] = val
end




"""
    Base.similar(V::T) where {T<:AbstractSpikeVector}

Create an (subtype of) [`LaskaCore.AbstractSpikeVector`](@ref) of the same type and size as `V` with the same samplerate.
"""
function Base.similar(V::T) where {T<:AbstractSpikeVector}
    T(similar(V.spiketimes), V.samplerate)
end

"""
    Base.similar(V::T, S::Type) where {T<:AbstractSpikeVector}

Create an (subtype of) [`LaskaCore.AbstractSpikeVector`](@ref) of the same size as `V` with type `S` with the same samplerate.
"""
function Base.similar(V::T, S::Type) where {T<:AbstractSpikeVector}
    T(similar(V.spiketimes, S), V.samplerate)
end



# Filter spiketimes

@inline Base.filter!(f::Function, V::SpikeVector) = Base.filter!(f, V.spiketimes)

@inline function Base.filter(f::Function, V::SpikeVector)
    out = deepcopy(V)
    Base.filter!(f, out)
    return out
end

# Unitful indexing and timerange filtering
# TODO: Implement/change these so '@views' can be used


function spikes_in_timerange(V::AbstractVector, lowerbound, upperbound)
    filter(inrange(lowerbound, upperbound), V)
end

function spikes_in_timerange(V::AbstractVector, range::AbstractRange)
    filter(inrange(range[begin], range[end]), V)
end

"""
    spikes_in_timerange(V::SpikeVector, range::AbstractRange{T}) where {T<:Quantity{<:Number,Unitful.𝐓}}

Returns a [`LaskaCore.SpikeVector`](@ref), only including spikes that fall within the times specified in `range`.
"""
function spikes_in_timerange(
    V::SpikeVector,
    range::AbstractRange{T},
) where {T<:Quantity{<:Number,Unitful.𝐓}}
    range = timetosamplerate(V, range)
    filter(inrange(ustrip(range[begin]), ustrip(range[end])), V)
end

"""
    spikes_in_timerange(V::SpikeVector, lowerbound::T, upperbound::T) where {T<:Quantity{<:Number,Unitful.𝐓}}

Returns a [`LaskaCore.SpikeVector`](@ref), only including spikes that fall within the times specified by `lowerbound` and `upperbound`.
"""
function spikes_in_timerange(
    V::SpikeVector,
    lowerbound::T,
    upperbound::T,
) where {T<:Quantity{<:Number,Unitful.𝐓}}
    lower = timetosamplerate(V, lowerbound)
    upper = timetosamplerate(V, upperbound)
    filter(inrange(lower, upper), V)
end

# In place versions

"""
    spikes_in_timerange!(V::AbstractVector, lowerbound, upperbound)

Filter a [`LaskaCore.SpikeVector`](@ref) in-place, only including spikes between `lowerbound` and `upperbound`.
"""
function spikes_in_timerange!(V::AbstractVector, lowerbound, upperbound)
    filter!(inrange(lowerbound, upperbound), V)
end

"""
    spikes_in_timerange!(V::SpikeVector, range::AbstractRange{T}) where {T<:Quantity{<:Number,Unitful.𝐓}}

Filter a [`LaskaCore.SpikeVector`](@ref) in-place using a [`Unitful`] time range, only including spikes that fall within the times specified in `range`.
"""
function spikes_in_timerange!(
    V::SpikeVector,
    range::AbstractRange{T},
) where {T<:Quantity{<:Number,Unitful.𝐓}}
    range = timetosamplerate(V, range)
    filter!(inrange(ustrip(range[begin]), ustrip(range[end])), V)
end

"""
    spikes_in_timerange!(V::SpikeVector, range::AbstractRange)

Filter a [`LaskaCore.SpikeVector`](@ref) in-place, only including spikes that fall within the times specified in `range`.
"""
function spikes_in_timerange!(V::SpikeVector, range::AbstractRange)
    filter!(inrange(range[begin], range[end]), V)
end

"""
    spikes_in_timerange!(V::SpikeVector, lowerbound::T, upperbound::T) where {T<:Quantity{<:Number,Unitful.𝐓}}

Filter a [`LaskaCore.SpikeVector`](@ref) in-place, only including spikes between `lowerbound` and `upperbound`.
"""
function spikes_in_timerange!(
    V::SpikeVector,
    lowerbound::T,
    upperbound::T,
) where {T<:Quantity{<:Number,Unitful.𝐓}}
    lower = timetosamplerate(V, lowerbound)
    upper = timetosamplerate(V, upperbound)
    filter!(inrange(lower, upper), V)
end

@inline inrange(lowerbound, upperbound) = x -> lowerbound <= x <= upperbound

#===============================
Functions on RelativeSpikeVector
===============================#

"""
    spikes_in_timerange(V::RelativeSpikeVector, lowerbound, upperbound)

Filter a [`LaskaCore.RelativeSpikeVector`](@ref), only including spikes between `lowerbound` and `upperbound`.
"""
function spikes_in_timerange(V::RelativeSpikeVector, lowerbound, upperbound)
    out = deepcopy(V)
    spikes_in_timerange!(out, lowerbound, upperbound)
    return out
end

"""
    spikes_in_timerange(V::RelativeSpikeVector, range::AbstractRange{T}) where {T<:Quantity{<:Number,Unitful.𝐓}}

Filter a [`LaskaCore.RelativeSpikeVector`](@ref) using a [`Unitful`] time range, only including spikes that fall within the times specified in `range`.
"""
function spikes_in_timerange(
    V::RelativeSpikeVector,
    range::AbstractRange{T},
) where {T<:Quantity{<:Number,Unitful.𝐓}}
    out = deepcopy(V)
    spikes_in_timerange!(out, range)
    return out
end

"""
    spikes_in_timerange(V::SpikeVector, lowerbound::T, upperbound::T) where {T<:Quantity{<:Number,Unitful.𝐓}}

Filter a [`LaskaCore.RelativeSpikeVector`](@ref) using [`Unitful`] bounds, only including spikes between `lowerbound` and `upperbound`.
"""
function spikes_in_timerange(
    V::RelativeSpikeVector,
    lowerbound::T,
    upperbound::T,
) where {T<:Quantity{<:Number,Unitful.𝐓}}
    out = deepcopy(V)
    spikes_in_timerange!(out, lowerbound, upperbound)
    return out
end

"""
    spikes_in_timerange(V::RelativeSpikeVector, range::AbstractRange)

Filter a [`LaskaCore.RelativeSpikeVector`](@ref), only including spikes that fall within the times specified in `range`.
"""
function spikes_in_timerange(V::RelativeSpikeVector, range::AbstractRange)
    out = deepcopy(V)
    spikes_in_timerange!(out, range)
    return out
end

# In-place versions for RelativeSpikeVector

"""
    spikes_in_timerange!(V::RelativeSpikeVector, lowerbound, upperbound)

Filter a [`LaskaCore.RelativeSpikeVector`](@ref) in-place, only including spikes between `lowerbound` and `upperbound`.
"""
function spikes_in_timerange!(V::RelativeSpikeVector, lowerbound, upperbound)
    for i in eachindex(V)
        spikes_in_timerange!(V[i], lowerbound, upperbound)
    end
end


"""
    spikes_in_timerange!(V::RelativeSpikeVector{T,U}, range::AbstractRange{TUNIT}) where {T,U,TUNIT<:Quantity{<:Number,Unitful.𝐓}}

Filter a [`LaskaCore.RelativeSpikeVector`](@ref) in-place using a [`Unitful`] time range, only including spikes that fall within the times specified in `range`.
"""
function spikes_in_timerange!(
    V::RelativeSpikeVector{T,U},
    range::AbstractRange{TimeUnit},
) where {T,U,TimeUnit<:Quantity{<:Number,Unitful.𝐓}}
    range = timetosamplerate(V, range)
    for i in eachindex(V)
        filter!(inrange(ustrip(range[begin]), ustrip(range[end])), V[i])
    end
end

"""
    spikes_in_timerange!(V::RelativeSpikeVector, range::AbstractRange)

Filter a [`LaskaCore.RelativeSpikeVector`](@ref) in-place, only including spikes that fall within the times specified in `range`.
"""
function spikes_in_timerange!(V::RelativeSpikeVector, range::AbstractRange)
    for i in eachindex(V)
        filter!(inrange(range[begin], range[end]), V[i])
    end
end

"""
    spikes_in_timerange!(V::RelativeSpikeVector, lowerbound::T, upperbound::T) where {T<:Quantity{<:Number,Unitful.𝐓}}

Filter a [`LaskaCore.RelativeSpikeVector`](@ref) in-place using [`Unitful`] bounds, only including spikes between `lowerbound` and `upperbound`.
"""
function spikes_in_timerange!(
    V::RelativeSpikeVector,
    lowerbound::T,
    upperbound::T,
) where {T<:Quantity{<:Number,Unitful.𝐓}}
    lower = timetosamplerate(V, lowerbound)
    upper = timetosamplerate(V, upperbound)
    for i in eachindex(V)
        filter!(inrange(lower, upper), V[i])
    end
end

# Generel RelativeSpikeVector interface functions

@inline function Base.setindex!(V::AbstractSpikeVector{T}, val::Vector{T}, i::Int) where {T}
    V.spiketimes[i] = val
end

"""
    (::Type{RelativeSpikeVector{T,U}})(::UndefInitializer, n, samplerate::U) where {T,U}

Create a [`LaskaCore.RelativeSpikeVector`](@ref) of length `n` with uninitialized sub-vectors.
"""
@inline function (::Type{RelativeSpikeVector{T,U}})(
    ::UndefInitializer,
    n::Int,
    sample_rate::U,
) where {T,U}
    RelativeSpikeVector{T,U}(Vector{Vector{T}}(undef, n), sample_rate)
end

# Creates a RelativeSpikeVector with each subvector[i] length of ns[i]
"""
    (::Type{RelativeSpikeVector{T,U}})(::UndefInitializer, ns::AbstractVector{Int}, sample_rate::U) where {T,U}

Create a [`LaskaCore.RelativeSpikeVector`](@ref) with `length(ns)` sub-vectors. Sub-vector`[n]` is initialized with `ns[n]` undefined elements. (ie `Vector{T}(undef, ns[n])`)
"""
@inline function (::Type{RelativeSpikeVector{T,U}})(
    ::UndefInitializer,
    ns::AbstractVector{Int},
    sample_rate::U,
) where {T,U}
    out = RelativeSpikeVector{T,U}(undef, length(ns), sample_rate)
    @inbounds for n in eachindex(out)
        out[n] = Vector{T}(undef, ns[n])
    end
    return out
end


@inline function (::Type{RelativeSpikeVector{T}})(
    ::UndefInitializer,
    n::Int,
    sample_rate::U,
) where {T,U}
    RelativeSpikeVector{T,U}(Vector{Vector{T}}(undef, n), sample_rate)
end

"""
    Base.similar(V::RelativeSpikeVector)

Create a [`LaskaCore.RelativeSpikeVector`](@ref) of the same type and dimensions as `V` and the same samplerate.
"""
function Base.similar(V::RelativeSpikeVector)
    RelativeSpikeVector(similar.(V), samplerate(V))
end

"""
    Base.similar(V::RelativeSpikeVector, S::Type)

Create a [`LaskaCore.RelativeSpikeVector`](@ref) of the same dimensions as `V` with type `S` and the same samplerate.
"""
function Base.similar(V::RelativeSpikeVector, S::Type)
    RelativeSpikeVector(similar.(V, S), samplerate(V))
end
