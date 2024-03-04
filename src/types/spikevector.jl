#=========================================================================
Vector type for holding spiketimes and the Hz/unit of their representation
=========================================================================#

struct SpikeVector{T,U} <: AbstractArray{T,1}
    spiketimes::Array{T,1}
    samplerate::U
end

# Interface implementation
Base.size(V::SpikeVector) = Base.size(V.spiketimes)

Base.getindex(V::SpikeVector, i::Int) = V.spiketimes[i]

Base.IndexStyle(::Type{<:SpikeVector}) = IndexLinear()

function Base.setindex!(V::SpikeVector{T}, v::T, i::Int) where {T}
    V.spiketimes[i] = v
end

Base.length(V::SpikeVector) = length(V.spiketimes)

Base.similar(V::SpikeVector) = SpikeVector(similar(V.spiketimes), V.samplerate)

samplerate(V::SpikeVector) = V.samplerate
