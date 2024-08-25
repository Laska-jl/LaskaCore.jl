# Convert a PhyOutput to a relativeSpikes struct


"""
    relativespikes(p::PhyOutput, stimtrain::Dict{String,T}, back, forward) where {T<:Real}

Returns a `RelativeSpikes` struct which contains `RelativeCluster`s and contains only spiketimes occuring `back` before triggers or `forward` after them.
Spiketimes are represented relative to the trigger event. For example a spike occuring 10ms before the trigger will appear as -300 if the sample rate
is 30 000Hz

# Arguments

All arguments representing time may be provided in a `Unitful.jl` unit or the samplerate of the experiment.
If one time is provided in a `Unitful` unit of time the others will have to be provided as such as well.

- `p::PhyOutput` -- A PhyOutput struct.
- `stimtrain` -- A `Dict` describing each trigger event. For example `Dict("US" => 0u"ms", "CS" => 300u"ms")` would mean one event (CS) at trigger and another (US) 300ms after trigger.
- `back` -- Time before trigger to include in each spike train.
- `forward` -- Time after trigger to include in each spike train.

"""
function relativespikes(p::PhyOutput{T,U,M,S}, stimtrain::Dict{String,V}, back::Y, forward::Y) where {T,U,V,M,S<:Vector{<:Integer},Y<:Real}
    # Create specs dict
    specs = (back=back, forward=forward, ntrig=ntrigs(p))
    # Convert back/forward to sample rate

    idvec = Vector{Int}(undef, 0)
    clustervec = Vector{RelativeCluster{T,U}}(undef, 0)
    for cluster in clustervector(p)
        samprate = samplerate(cluster)
        resvec = RelativeSpikeVector{T}(undef, length(triggertimes(p)), samplerate(cluster))
        for n in eachindex(resvec)
            resvec[n] = T[]
        end
        _filtertriggers!(resvec, spiketimes(cluster), triggertimes(p), back, forward)
        push!(
            clustervec,
            RelativeCluster(id(cluster), info(cluster), resvec)
        )
        push!(idvec, id(cluster))
    end
    return RelativeSpikes(idvec, clustervec, triggertimes(p), getmeta(p), info(p), stimtrain, specs)
end

function relativespikes(p::PhyOutput, stimtrain::Dict, back::TUnit, forward::TUnit)
    backsamp = timetosamplerate(p, back)
    forwardsamp = timetosamplerate(p, forward)
    stimtrainsamp = Dict(key => timetosamplerate(p, val) for (key, val) in pairs(stimtrain))
    relativespikes(p, stimtrainsamp, backsamp, forwardsamp)
end

function _filtertriggers!(resultvector::RelativeSpikeVector{T,U}, spiketimes::AbstractVector{T}, triggertimes, back, forward) where {T<:Real,U}
    pos = 1
    minT = triggertimes[pos] - back
    maxT = triggertimes[pos] + forward
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
