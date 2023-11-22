

# Extract triggers from Vector of triggerchannel


# TODO: Update this and put in readglx.jl
function getchan(
    p::PhyOutput,
    ch::Union{Int,Vector{Int},UnitRange{Int64}},
    tmin::Union{Float64,Int},
    tmax::Union{Float64,Int,String},
    converttoV::Bool=true,
    threading::Bool=false
)

    filemax::Int64 = Int64(parse(Int64, p._meta["fileSizeBytes"]) / parse(Int64, p._meta["nSavedChans"]) / 2)
    samplfrq::Float64 = parse(Float64, p._meta["imSampRate"])
    # Convert times (IN SECONDS) to sample frequencies
    # Add 1 since indexing is 1-based(?)

    #Check if tmin < 0
    if tmin < 0
        ArgumentError("tmin must be above 0")
    end

    # Check if max time should be used
    if typeof(tmax) == String
        if tmax != "max"
            ArgumentError("tmax must be a number or the string 'max'")
        end
        tmax = filemax
    else
        # Convert tmax to sample frequency
        tmax::Int64 = Int(round((tmax * samplfrq))) + 1::Int64
    end
    # Check if tmin >= tmax
    if tmin >= tmax
        ArgumentError("tmin must be <tmax")
    end

    # Convert tmin to sample frequency
    tmin::Int64 = Int(round((tmin * samplfrq))) + 1::Int64

    if tmax > filemax
        ArgumentError("tmax larger than length of recording")
    end

    # Memory map data and load the selected chunks
    karta::Matrix{Int16} = spikemmap(p)
    len::UnitRange{Int64} = tmin:tmax
    #it::Vector{Tuple{Int64, Int64}} = collect(enumerate(len))

    if threading
        r::Union{Matrix{Int16},Vector{Int16}} = Matrix{Int16}(undef, length(len), length(ch))
        t::Vector{UnitRange{Int64}} = collect(Iterators.partition(len, 10000))
        n::Vector{UnitRange{Int64}} = collect(Iterators.partition(1:length(len), 10000))
        it = tuple.(n, t)
        #it = collect(enumerate(len))


        if length(ch) > 1
            Laska.importchx!(ch, r, karta, it)
        elseif length(ch) == 1
            Laska.importch1!(ch, r, karta, it)
        end
    else
        r = karta[ch, len]
    end

    if converttoV
        conv::Union{Matrix{Float64},Vector{Float64}} = tovolts(p._meta, r)
        return conv
    else
        return r
    end

end # Getchan



function importchx!(channels::Union{Int,Vector{Int},UnitRange{Int64}}, a::Matrix{Int16}, mm::Matrix{Int16}, i::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}})
    for (ntim, t) in i
        a[ntim, :] = transpose(mm[channels, t])
    end
    return a
end

function importch1!(channels::Union{Int,Vector{Int},UnitRange{Int64}}, a::Vector{Int16}, mm::Matrix{Int16}, i::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}})
    for (ntim, t) in i
        a[ntim, 1] = mm[channels, t]
    end
    return a
end

function spikemmap(p::PhyOutput)
    n::Int = parse(Int, getmeta(p, "nSavedChans"))
    s::Int = Int(parse(Int, getmeta(p, "fileSizeBytes")) / (2 * n))
    tmp::IOStream = open(p._binpath, "r")
    m::Array{Int16,2} = mp.mmap(tmp, Array{Int16,2}, (n, s), 0)
    close(tmp)
    return m
end # spikemmap

function tovolts(meta::Dict{SubString{String},SubString{String}}, i::Union{Vector{Int16},Array{Int16}})
    Imax::Float64 = parse(Float64, meta["imMaxInt"])
    Vmax::Float64 = parse(Float64, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl = meta["~imroTbl"]
        tbl = split(tbl, ")(")
        tbl = split(tbl[2], " ")
        gain = parse(Int, tbl[4])
    else
        gain::Float64 = 80
    end
    cfactor::Float64 = Vmax / Imax / gain
    return i .* cfactor
end # tovolts

function importchanint16(path::String="")
    if path == ""
        path = Gtk.open_dialog_native("Select trigger file", action=GtkFileChooserAction.GTK_FILE_CHOOSER_ACTION_OPEN)
    end
    if path[end-3:end] == ".bin"
        tmp::IOStream = open(path, "r")
        res = Vector(reinterpret(Int16, read(tmp)))
        close(tmp)
        return res
    elseif path[end-3:end] == ".csv"
        res2::Vector = Vector(Matrix(CSV.read(path, DataFrame, header=false))[:, 1])
        return res2
    else
        ArgumentError("Only '.bin' or '.csv' files allowed.")
    end
end #importchan
