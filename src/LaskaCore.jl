module LaskaCore

using NPZ
using CSV
using DataFrames
using Mmap

include("types/cluster.jl")
# Export AbstracCluster and related methods
export AbstractCluster, id, nspikes, info, spiketimes
# Export Cluster and related types
export Cluster
export RelativeCluster
include("types/experiments.jl")
# Export phyoutput
export AbstractExperiment, getcluster, ntrigs, clusterids, triggertimes, clustervector, getmeta
export PhyOutput
export RelativeSpikes, relativespecs, stimtimes

# Importing
include("import/importphy.jl")
# Export function for importing phy output
export importphy
include("import/importmisc.jl")
include("import/readglx.jl")
# Export glx reading functions
export spikemmap, tovolts, tovolts!, parseglxmeta

end
