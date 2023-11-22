#---------------------------------------------#
# Functions for importing spikeGLX .bin files |
#---------------------------------------------#

"""

    spikemmap(file::String, nchans::Int, filesizebytes::Int)
    spikemmap(file::String, meta::Dict{Substring{String}, SubString{String}})

Create a memory map of a spikeGLX .bin file. Requires a path to the `file`, the number of channels and the file size in bytes.              
The easiest way to provide this is to pass a parsed .meta file.

"""
function spikemmap(file::String, nchans::Int, filesizebytes::Int)
    tmp::IOStream = open(file, "r")
    mm::Array{Int16,2} = Mmap.mmap(tmp, Array{Int16,2}, (nchans, filesizebytes), 0)
    close(tmp)
    return mm
end

function spikemmap(file::String, meta::Dict{SubString{String},SubString{String}})
    n::Int = parse(Int, meta["nSavedChans"])
    s::Int = Int(parse(Int, meta["fileSizeBytes"]) / (2 * n))
    return spikemmap(file, n, s)
end


"""
    tovolts(in::Matrix{Int16}, meta::Dict{SubString{String},SubString{String}})

Convert a Vector/Matrix of raw spikeGLX `Int16` data to volts.

"""
function tovolts(in::Matrix{Int16}, meta::Dict{SubString{String},SubString{String}})
    out::Matrix{Float32} = similar(in)
    return tovolts!(out, in, meta)
end

function tovolts!(out::Matrix{Float32}, in::Matrix{Int16}, meta::Dict{SubString{String},SubString{String}})
    Imax::Float32 = parse(Float32, meta["imMaxInt"])
    Vmax::Float32 = parse(Float32, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl = meta["~imroTbl"]
        tbl = split(tbl, ")(")
        tbl = split(tbl[2], " ")
        gain::Float32 = parse(Float32, tbl[4])
    else
        gain = 80.0
    end
    cfactor::Float32 = Vmax / Imax / gain
    return tovolts!(out, in, cfactor)
end


function tovolts!(out::Matrix{Float32}, in::Matrix{Int16}, cfactor::Float32)
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
    res::Vector{Int16} = Vector(reinterpret(Int16, read(tmp)))
    close(tmp)
    return res
end

"""

    importchanint16csv(path::String)

Import a .csv file containing a **SINGLE** channel.

"""
function importchanint16csv(path::String)
    res::Vector{Int16} = Vector(Matrix(CSV.read(path, DataFrame, header=false))[:, 1])
    return res

end

"""

    gettrig(t::Vector{T}) where {T<:Real}

Extract trigger event indices.

Assumes that `t` is 0 except at trigger events. Returns only the *first* index at each trigger event(ie series of sequential nonzero indices).

"""
function gettrig(t::Vector{T}) where {T<:Real}
    r::Vector{Int64} = findall(!iszero, t)
    p::Matrix{Int64} = hcat(getindex.(r, 1), getindex.(r - circshift(r, 1), 1))
    s::Vector{Int64} = @inbounds p[p[:, 2].!=1, 1]
    return s
end

"""

    parseglxmeta(file::String)

Parses a spikeGLX .meta file into a Dict.

"""
function parseglxmeta(file::String)
    tmp = open(file, "r")
    metaraw = readlines(tmp)
    close(tmp)
    metaraw = split.(metaraw, "=")
    return Dict{SubString{String},SubString{String}}(i[1] => i[2] for i in metaraw)
end
