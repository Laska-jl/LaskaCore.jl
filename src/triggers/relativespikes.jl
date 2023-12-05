################################################
# 
# Convert a PhyOutput to a relativeSpikes struct
#
################################################
# DONE: Kom på ett bättre namn


"""
    relativespikes(p::PhyOutput{T}, stimtrain::Dict{String,T}, back::T, forward::T) where {T<:Real}

Returns a `RelativeSpikes` struct which wraps `RelativeCluster`:s and contains only spiketimes occuring `back` ms before triggers or `forward` ms after them.

# Variables
- `p::PhyOutput{T}` -- A PhyOutput struct.
- `stimtrain::Dict{String, T}` -- A dict describing each trigger event. For example `Dict("CS" => 0, "US" => 300)` would mean one event(CS) at trigger and another(US) 300ms after trigger.
- `back::T` -- Number of ms before trigger to include.
- `forward::T` -- Number of ms after trigger to include.
"""
function relativespikes(p::PhyOutput{T}, stimtrain::Dict{String,T}, back::T, forward::T) where {T<:Real}
    # Create specs dict
    specs = Dict{String,T}(
        "ntrig" => ntrigs(p),
        "back" => back,
        "forward" => forward
    )
    # Convert back/forward to sample rate
    backF::T = mstosamplerate(back, p)
    forwardF::T = mstosamplerate(forward, p)

    idvec::Vector{Int} = Vector{Int}(undef, 0)
    clustervec::Vector{RelativeCluster{T}} = Vector{RelativeCluster{T}}(undef, 0)
    for cluster in clustervector(p)
        resvec = Vector{Vector{T}}(undef, length(triggertimes(p)))
        for n in eachindex(resvec)
            resvec[n] = T[]
        end
        filtertriggers!(resvec, spiketimes(cluster), triggertimes(p), backF, forwardF)
        push!(
            clustervec,
            RelativeCluster(id(cluster), info(cluster), resvec)
        )
        push!(idvec, id(cluster))
    end
    return RelativeSpikes(idvec, clustervec, triggertimes(p), getmeta(p), info(p), stimtrain, specs)
end


function filtertriggers!(resultvector::Vector{Vector{T}}, spiketimes::Vector{T}, triggertimes::Vector{T}, back::T, forward::T) where {T<:Real}
    pos::Int = 1
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
