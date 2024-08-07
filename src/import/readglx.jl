#---------------------------------------------#
# Functions for importing spikeGLX .bin files |
#---------------------------------------------#

# function importglx(path::String, meta::Dict, channels, t::U) where {U<:LaskaCore.TUnit}
#     nchans = meta["nSavedChans"]
#     
#
# end

const NP1_LIKE_PROBES::NTuple{10,Int64} = (0, 1020, 1030, 1100, 1120, 1121, 1122, 1123, 1200, 1300)
const NP2_SINGLE_MULTIPLEX_SHANK_PROBES::NTuple{3,Int64} = (21, 2003, 2004)
const NP2_4_SHANK_PROBES::NTuple{3,Int64} = (24, 2013, 2014)
const NP_TYPE_2020_QUAD_PROBES::NTuple{1,Int64} = (2020,)
const NP_UHD_PROGRAMMABLE_PROBES::NTuple{1,Int64} = (1110,)

abstract type AbstractImroTbl{ProbeType} end
abstract type AbstractNPProbe end

struct NP1Like <: AbstractNPProbe end

function probetypes(::Type{NP1Like})
    Set(NP1_LIKE_PROBES)
end

struct NP2Single <: AbstractNPProbe end

function probetypes(::Type{NP2Single})
    Set(NP2_SINGLE_MULTIPLEX_SHANK_PROBES)
end

struct NP24Shank <: AbstractNPProbe end

function probetypes(::Type{NP24Shank})
    Set(NP2_4_SHANK_PROBES)
end

struct NPQuad <: AbstractNPProbe end

function probetypes(::Type{NPQuad})
    Set(NP_TYPE_2020_QUAD_PROBES)
end

struct NPUHDProgrammable <: AbstractNPProbe end

function probetypes(::Type{NPUHDProgrammable})
    Set(NP_UHD_PROGRAMMABLE_PROBES)
end

"""
    struct ImroTbl{T<:Real,P<:AbstractNPProbe,E<:NamedTuple} <: AbstractImroTbl{T}
        table::Matrix{T}
        entries::E
    end

Struct returned by [`LaskaCore.parseimroTbl`](@ref). Contains a matrix of the entries in the
Imec Readout Table and a named tuple of the different entries. Each column is a different entry
(such as Channel ID, AP band gain etc) while each row contain the values of 1 channel.

For more information on possible entries depending on probetype please see the [SpikeGLX manual](http://billkarsh.github.io/SpikeGLX/Sgl_help/Metadata_30.html)

In order to get the entries of a certain parse `ImroTbl` `Base.keys()` may be used.

# Example

The table may be indexed to as such:
```julia
table = LaskaCore.parseimroTbl(SomePhyOutput)

# Will return the entire column of channel ids
table[:Channel_id]

# Will return the channel id of the first channel
table[1,:Channel_id]

# Will return entries for the first channel
table[1,:]

# Will return the 4th entry of the 1st channel
table[1, 4]

```
"""
struct ImroTbl{T<:Real,P<:AbstractNPProbe,E<:NamedTuple} <: AbstractImroTbl{T}
    table::Matrix{T}
    entries::E
end

function Base.keys(t::AbstractImroTbl)
    keys(t.entries)
end

"""
    parseimroTbl(tbl)

Determines the probe type and parses the Imec Readout Table accordingly.

For information on what info is contained in the table depending on probe type
please see the section on `~imroTbl` in the [SpikeGLX manual](http://billkarsh.github.io/SpikeGLX/Sgl_help/Metadata_30.html)

Currently, parsing functions for the following probes are implemented:          

- NP 1.0-like (0,1020,1030,1100,1120,1121,1122,1123,1200,1300)

Returns a [`LaskaCore.ImroTbl`](@ref)

"""
function parseimroTbl(tbl)
    s = split(tbl[begin+1:end-1], ")(")
    ptype, nchans = parse.(Int64, split(s[1], ","))
    if ptype in probetypes(NP1Like)
        __parseimroTbl(s, nchans, NP1Like)
    else
        throw(ErrorException("Parsing function for probe $ptype not implemented"))
    end
end

function parseimroTbl(p::PhyOutput)
    parseimroTbl(getmeta(p, "~imroTbl"))
end

function __parseimroTbl(s, nchans::Integer, ::Type{NP1Like})
    outm = zeros(Int64, nchans, 6)
    entries = (
        Channel_id=1,
        BankNr=2,
        RefIdInd=3,
        APGain=4,
        LFGain=5,
        APHiPass=6
    )
    for r in Iterators.drop(eachindex(s), 1)
        outm[r-1, :] = parse.(Int64, split(s[r, :][1], " "))
    end
    return ImroTbl{Int64,NP1Like,@NamedTuple{
        Channel_id::Int64,
        BankNr::Int64,
        RefIdInd::Int64,
        APGain::Int64,
        LFGain::Int64,
        APHiPass::Int64
    }}(outm, entries)
end

function Base.getindex(t::AbstractImroTbl, e::Symbol)
    t.table[:, t.entries[e]]
end

function Base.getindex(t::AbstractImroTbl, r::Integer, e::Symbol)
    t.table[r, t.entries[e]]
end

function Base.getindex(t::AbstractImroTbl, r, c)
    t.table[r, c]
end


"""
    importglx(path::String, channels, t, meta)

Import raw Neuropixels channels from a .ap.bin file.

# Arguments

- `path`: Direct path to the ap.bin file from which to read
- `channels`: Channels to import. May be any value capable of indexing rows of a matrix.
- `t`: Time indices to include. Using 1:30001 will import the full first second of the recording if the samplerate is 30000 Hz.
- `meta`: A parsed SpikeGLX .ap.meta file as imported using [`LaskaCore.importphy`](@ref) or [`LaskaCore.parseglxmeta`](@ref).

In order to convert the raw data to volts the function [`LaskaCore.tovolts`](@ref) can be used.

"""
function importglx(path::String, channels, t, meta)
    mm = spikemmap(path, meta)
    return mm[channels, t]
end

"""
    spikemmap(file::String, nchans::Int)
    spikemmap(file::String, meta::Dict{String, String})

Create a memory map of a spikeGLX .bin file. Requires a path to the `file` and the number of channels.              
The easiest way to provide the number of channels is to pass a parsed .meta file.

"""
function spikemmap(file::String, nchans::Int)
    filesizebytes = filesize(file)
    tmp::IOStream = open(file, "r")
    mm::Array{Int16,2} = Mmap.mmap(tmp, Array{Int16,2}, (nchans, Int(filesizebytes / (2 * nchans))), 0)
    close(tmp)
    return mm
end

function spikemmap(file::String, meta::Dict{String,String})
    n::Int = parse(Int, meta["nSavedChans"])
    return spikemmap(file, n)
end


"""
    tovolts(in::Matrix{Int16}, meta::Dict{String,String})

Convert a Vector/Matrix of a raw spikeGLX AP recording to volts.

The meta expected is the `Dict` returned from [`LaskaCore.getmeta`](@ref) or [`LaskaCore.parseglxmeta`](@ref).

"""
function tovolts(in::AbstractArray, meta::Dict{String,String})
    # TODO: Fix these so they can use vectors and arrays
    out = similar(in, Float64)
    return tovolts!(out, in, meta)
end

function tovolts!(out::AbstractArray, in::AbstractArray, meta::Dict{String,String})
    Imax = parse(Float64, meta["imMaxInt"])
    Vmax = parse(Float64, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl = parseimroTbl(meta["~imroTbl"])
        gains = unique(tbl[:APGain])
        length(gains) != 1 && throw(ErrorException("Multiple gains found in Imec Readout Table"))
        gain = Float64(gains[1])
    else
        gain = 80.0
    end
    cfactor = Vmax / Imax / gain
    return tovolts!(out, in, cfactor)
end


function tovolts!(out::AbstractArray{<:AbstractFloat}, in::AbstractArray, cfactor::AbstractFloat)
    if size(out) != size(in)
        throw(ArgumentError("out and in matrices should have the same dimensions!"))
    end
    for i in eachindex(out)
        out[i] = in[i] * cfactor
    end
    return out
end

"""
    importchanint16bin(path::String)

Import a .bin file containing a **SINGLE** channel.

"""
function importchanint16bin(path::String)
    tmp::IOStream = open(path, "r")
    res = Vector{Int16}(reinterpret(Int16, read(tmp)))
    close(tmp)
    return res
end

"""

    importchanint16csv(path::String)

Import a .csv file containing a **SINGLE** channel. This is considerably slower than importing
from a .bin-file as exported by SpikeGLX.

"""
function importchanint16csv(path::String)
    res::Vector{Int16} = Vector(Matrix(CSV.read(path, DataFrame, header=false))[:, 1])
    return res
end

"""

    gettrig(t::Vector{T}) where {T<:Real}

Extract trigger event indices.

Assumes that `t` is 0 except at trigger events. Returns only the *first* index at each trigger event(ie series of sequential nonzero indices).

# Example

```julia

# * = Indices that will be returned, ^ = trigger
#            *                 *
v = [0,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0]
#            ^ ^ ^ ^           ^ ^ ^ ^ ^ ^ 

# Will return [5, 14]
LaskaCore.gettrig(v)
```

"""
function gettrig(t::Vector{T}) where {T<:Real}
    # Find all indices with triggers
    r = findall(!iszero, t)
    # Prepare result buffer
    c = Vector{Int64}(undef, length(r))
    n = 1
    # Add the first index of r to c since it will always mark the start of an event.
    c[begin] = r[begin]
    # Iterate over all indices except first since it's already covered,
    @views for i in Iterators.drop(eachindex(r), 1)
        # If indices are nonsequential the current one must mark the start of an event.
        if r[i] - r[i-1] != 1
            n += 1
            c[n] = r[i]
        end
    end
    return c[begin:n]
end

#function gettrig(t::Vector{T}) where {T<:Real}
#    r = findall(!iszero, t)
#    p::Matrix{Int64} = hcat(getindex.(r, 1), getindex.(r - circshift(r, 1), 1))
#    s::Vector{Int64} = @inbounds p[p[:, 2].!=1, 1]
#    return s
#end

"""

    parseglxmeta(file::String)

Parses a spikeGLX .meta file into a Dict.

"""
function parseglxmeta(file::String)
    tmp = open(file, "r")
    metaraw = readlines(tmp)
    close(tmp)
    metaraw = split.(metaraw, "=")
    return Dict{String,String}(String(i[1]) => String(i[2]) for i in metaraw)
end

