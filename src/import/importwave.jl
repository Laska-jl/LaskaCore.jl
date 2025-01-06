function importwaves(binp::String, sp::AbstractVector, meta::Dict{String,String}, win::Tuple{<:Int,<:Int})
    mm = spikemmap(binp, meta)

    Imax = parse(Float32, meta["imMaxInt"])
    Vmax = parse(Float32, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl = parseimroTbl(meta["~imroTbl"])
        gains = unique(tbl[:APGain])
        length(gains) != 1 && throw(ErrorException("Multiple gains found in Imec Readout Table"))
        gain = Float32(gains[1])
    else
        gain = 80.0f0
    end
    cfactor = Vmax / Imax / gain

    n_AP_chans = parse(Int64, parsesnsChanMap(meta["~snsChanMap"])[1,1])

    out = Array{Float32}(undef, n_AP_chans, win[1] + win[2] + 1, length(sp))

    for i in eachindex(sp)
        out[:, :, i] = mm[begin:n_AP_chans, sp[i]-win[1]:sp[i]+win[2]]
    end
    for i in axes(out, 3)
        tovolts!(@view(out[:, :, i]), @view(out[:, :, i]), cfactor)
    end

    return out
end


