# Convert a PhyOutput to a relativeSpikes struct


"""
    relativespikes(p::PhyOutput{T}, stimtrain::Dict{String,T}, back, forward) where {T<:Real}

Returns a `RelativeSpikes` struct which wraps `RelativeCluster`:s and contains only spiketimes occuring `back` ms before triggers or `forward` ms after them.

# Variables
- `p::PhyOutput{T}` -- A PhyOutput struct.
- `stimtrain::Dict{String, T}` -- A dict describing each trigger event. For example `Dict("CS" => 0, "US" => 300)` would mean one event(CS) at trigger and another(US) 300ms after trigger.
- `back::T` -- Number of ms before trigger to include in each spike train.
- `forward::T` -- Number of ms after trigger to include in each spike train.
"""
function relativespikes(p::PhyOutput{T,U,M,S}, stimtrain::Dict{String,V}, back::Y, forward::Y) where {T,U,V,M,S<:Vector{<:Integer},Y<:Real}
    # Create specs dict
    specs = (back=back, forward=forward, ntrig=ntrigs(p))
    # Convert back/forward to sample rate

    idvec = Vector{Int}(undef, 0)
    clustervec = Vector{RelativeCluster{T,U}}(undef, 0)
    for cluster in clustervector(p)
        samprate = samplerate(cluster)
        backF = back * samprate / 1000
        forwardF = forward * samprate / 1000
        resvec = RelativeSpikeVector{T}(undef, length(triggertimes(p)), samplerate(cluster))
        for n in eachindex(resvec)
            resvec[n] = T[]
        end
        _filtertriggers!(resvec, spiketimes(cluster), triggertimes(p), backF, forwardF)
        push!(
            clustervec,
            RelativeCluster(id(cluster), info(cluster), resvec)
        )
        push!(idvec, id(cluster))
    end
    return RelativeSpikes(idvec, clustervec, triggertimes(p), getmeta(p), info(p), stimtrain, specs)
end


function _filtertriggers!(resultvector::RelativeSpikeVector{T,U}, spiketimes::AbstractVector{T}, triggertimes, back, forward) where {T<:Real,U}
    pos = 1
    minT::T = triggertimes[pos] - back
    maxT::T = triggertimes[pos] + forward
    for n in eachindex(spiketimes)
        if minT <= spiketimes[n] <= maxT
            push!(resultvector[pos], spiketimes[n] - triggertimes[pos])
        elseif spiketimes[n] > maxT
            pos += 1
            if pos > length(triggertimes)
                break
            end
            minT = triggertimes[pos] - back
            maxT = triggertimes[pos] + forward
        end
    end
end


function relativespikes(p::PhyOutput, stimtrain::Dict, back::T, forward::T) where {T<:LaskaCore.TUnit}
    backms = Unitful.uconvert(u"ms", back) |> ustrip
    forwardms = Unitful.uconvert(u"ms", forward) |> ustrip

    relativespikes(p, stimtrain, backms, forwardms)
end
