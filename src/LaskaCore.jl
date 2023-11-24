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

# Helper functions
include("helpers/timeconv.jl")
include("helpers/rounding.jl")
include("helpers/spikesatdepth.jl")
export spikesatdepth
include("helpers/normalize.jl")
include("helpers/unpackvector.jl")
include("helpers/isntempty.jl")
include("helpers/findmax.jl")
include("helpers/isi.jl")

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

# Summarizing statistics
include("summarize/cv2.jl")
include("summarize/mad.jl")
include("summarize/frequency.jl")
include("summarize/relativefrequency.jl")
export cv2, cv2mean, mad, relativefrequency, frequency

end
