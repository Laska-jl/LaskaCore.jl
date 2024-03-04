module LaskaCore

using NPZ
using CSV
using DataFrames
using Mmap
using Dates
using Unitful

include("types/cluster.jl")
# Export AbstracCluster and related methods
export AbstractCluster, id, nspikes, info, spiketimes, timetosamplerate
# Export Cluster and related types
export Cluster
export RelativeCluster
include("types/experiments.jl")
# Export phyoutput
export AbstractExperiment,
    getcluster, ntrigs, clusterids, triggertimes, clustervector, getmeta
export PhyOutput
export RelativeSpikes, relativespecs, stimtimes

include("types/spikevector.jl")
export SpikeVector

# Helper functions
include("helpers/timeconv.jl")
include("helpers/rounding.jl")
include("helpers/spikesatdepth.jl")
export spikesatdepth
include("helpers/unpackvector.jl")
include("helpers/isntempty.jl")
include("helpers/findmax.jl")

# Work around triggers
include("triggers/relativespikes.jl")
export relativespikes

# Importing
include("import/importphy.jl")
# Export function for importing phy output
export importphy
include("import/importmisc.jl")
include("import/readglx.jl")
# Export glx reading functions
export spikemmap, tovolts, tovolts!, parseglxmeta

include("filter/filterexperiment.jl")


end
