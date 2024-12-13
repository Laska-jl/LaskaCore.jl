


"""
    importwaves(binp, sp, meta, win[, ch])

Import raw waveforms surrounding spiketimes `sp` from a `.ap.bin` file.

# Arguments

- `binp`: Path to the `.ap.meta` file.
- `sp`: Vector of spike _indices_ of interest. A spike at 1.0s in a file recorded at 30 000 Hz would have an index of 30000.
- `meta`: Meta dictionary built from the appropriate `.ap.meta` file. Can be imported along with with [`LaskaCore.importphy`](@ref) or [`LaskaCore.parsglxmeta`](@ref).
- `win`: A `Tuple{T,T}` specifying the length of the window (in terms of indices) around each spike in `sp` to include.
- `ch`: Optional argument specifying channel(s) to include. If omitted all AP channels are included. May be an `Integer`, `AbstractUnitRange` or `AbstractVector{<:Integer}`.
"""
function importwaves end

function importwaves(binp::String, sp::AbstractVector{T}, meta::Dict{String,String}, win::Tuple{<:T,<:T}) where {T<:Integer}
    mm = spikemmap(binp, meta)

    Imax = parse(Float32, meta["imMaxInt"])
    Vmax = parse(Float32, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl = parseimroTbl(meta["~imroTbl"])
        gains = unique(tbl[:APGain])
        length(gains) != 1 && throw(ErrorException("Multiple gains found in Imec Readout Table, this is not supported at this time"))
        gain = Float32(gains[1])
    else
        gain = 80.0f0
    end
    cfactor = Vmax / Imax / gain

    n_AP_chans = parse(Int64, parsesnsChanMap(meta["~snsChanMap"])[1, 1])

    # out = Array{Float32}(undef, n_AP_chans, win[1] + win[2] + 1, length(sp))

    # for i in eachindex(sp)
    inds = __solve_indices(sp, win)
    raw = mm[begin:n_AP_chans, inds]

    out = Float32.(reshape(raw, (n_AP_chans, sum(win) + 1, length(sp))))

    for i in axes(out, 3)
        tovolts!(@view(out[:, :, i]), @view(out[:, :, i]), cfactor)
    end

    return out
end


function importwaves(binp::String, sp::AbstractVector{T}, meta::Dict{String,String}, win::Tuple{<:T,<:T}, ch::Union{<:Integer,<:AbstractUnitRange,<:AbstractVector{<:Integer}}) where {T<:Integer}
    mm = spikemmap(binp, meta)

    Imax = parse(Float32, meta["imMaxInt"])
    Vmax = parse(Float32, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl = parseimroTbl(meta["~imroTbl"])
        gains = unique(tbl[:APGain])
        length(gains) != 1 && throw(ErrorException("Multiple gains found in Imec Readout Table, this is not supported at this time"))
        gain = Float32(gains[1])
    else
        gain = 80.0f0
    end
    cfactor = Vmax / Imax / gain

    n_AP_chans = parse(Int64, parsesnsChanMap(meta["~snsChanMap"])[1, 1])

    inds = __solve_indices(sp, win)
    raw = mm[ch, inds]

    out = Float32.(reshape(raw, (length(ch), sum(win) + 1, length(sp))))

    for i in axes(out, 3)
        tovolts!(@view(out[:, :, i]), @view(out[:, :, i]), cfactor)
    end

    return out
end


# Returns a vector of indices which, when applied to an array/vector
# returns the data sp - win[1]:sp + win[1] for each element of sp
function __solve_indices(sp::AbstractVector{T}, win::Tuple{T,T}) where {T<:Integer}
    winlen = sum(win) + 1
    out = Vector{T}(undef, length(sp) * winlen)
    idx_start = 1
    for i in eachindex(sp)
        out[idx_start:idx_start+winlen-1] .= sp[i]-win[1]:sp[i]+win[2]
        idx_start += winlen
    end
    return out
end
