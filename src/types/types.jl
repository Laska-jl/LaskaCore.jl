struct ClusterBaseline
    x::Dict{Int64,Matrix{Float64}}
    function ClusterBaseline(c::Dict{Int64,Matrix{Float64}})
        new(c)
    end
end

function Base.keys(clusterbaseline::Laska.ClusterBaseline)
    return keys(clusterbaseline.x)
end

Base.show(io::IO, clusterbaseline::ClusterBaseline) = println(
    "Pre stimulus baselines by cluster\nClusters: $(keys(clusterbaseline))"
)

struct DepthBaseline
    x::Dict{Int64,Matrix{Float64}}
    function DepthBaseline(d::Dict{Int64,Matrix{Float64}})
        return new(d)
    end
end

function Base.keys(depthbaseline::Laska.DepthBaseline)
    return keys(depthbaseline.x)
end

Base.show(io::IO, depthbaseline::DepthBaseline) = println(
    "Pre stimulus baselines by depth\nDepths: $(keys(depthbaseline))"
)



