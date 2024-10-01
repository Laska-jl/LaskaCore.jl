

function extractwaves(raw::Matrix{T}, sp::SpikeVector, raw_start::Integer, win_b::Integer, win_f::Integer) where {T}
    p = spikes_in_timerange(sp, raw_start, raw_start + size(raw, 2))
    out = Array{T}(undef, size(raw, 1), win_b + win_f + 1, length(p))
    for i in eachindex(p)
        out[:, :, i] = raw[:, p[i]-win_b+raw_start:p[i]+win_f+raw_start]
    end
    return out
end

function filterwaves(waves::Array{T,3}, responsetype::FilterType=Highpass(300, fs=30000), filterdesign::FilterCoefficients=Butterworth(3)) where {T}
    out = similar(waves)
    coefs = digitalfilter(responsetype, filterdesign)
    Threads.@threads for k in axes(waves, 3)
        for i in axes(waves, 1)
            out[i, :, k] = filtfilt(coefs, waves[i, :, k])
        end
    end
    return out
end

function avgwaves(waves::Array{T,3}, n::Int = 10) where {T}
    chs, samples, n_waves = size(waves)
    new = div(n_waves, n) + (mod(n_waves, n) == 0 ? 0 : 1)
    out = Array{T,3}(undef, chs, samples, new)
    for i in 1:n:n_waves
        out[:,:,div(i, n)+1] = mean(waves[:,:,i:min(i+n-1, n_waves)], dims=3)
    end
    return out
end

# https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7951023

function alignwave(wave::AbstractVector, len::Int)
    out = similar(wave)
    out[begin] = wave[begin]
    for n in Iterators.take(eachindex(wave), length(wave)-1)
        L = min(n, len)
        out[n + 1] = out[n] + __o(wave, n, L) +
            (-2 / L) * (wave[n+1] - wave[n - L + 1]) + wave[n + 1] - wave[n - L + 1]
    end
    return out
end

function __o(in, n, L)
    -2/L * sum([in[n - i] for i in 0:L-1])
end


function alignwave2(wave::AbstractVector, len::Int)
    out = deepcopy(wave)
    L=len
    for n in Iterators.drop(Iterators.take(eachindex(wave), length(wave)-1), len + 1)
        out[n + 1] = out[n] + __o2(wave, n, L) +
            (-2 / L) * (wave[n+1] - wave[n - L - 1]) + wave[n + 1] - wave[n - L]
    end
    return out
end

function __o2(in, n, L)
    -2/L * sum([in[n - i] for i in 0:L-1])
end

function wave_baseline(wave::AbstractVector, n_samples::Int)
    m = (mean(@view(wave[begin:begin+n_samples])) + mean(@view(wave[end-n_samples:end]))) / 2
end

function adjust_wave_baseline!(wave::AbstractVector, n_samples::Int = 10)
    wave .-= wave_baseline(wave, n_samples)
end

"""
    adjust_wave_baselines!(waves::AbstractMatrix{T}, n_samples::Int=10) where {T}

Adjusts the baseline of a matrix (n_channels x n_samples) of waveforms by subtracting the 
mean of the first and last `n_samples` from each row. Each channel is adjusted individually.

"""
function adjust_wave_baselines!(waves::AbstractMatrix{T}, n_samples::Int=10) where {T}
    ms = (mean(@view(waves[:, begin:begin+n_samples]), dims=2) .+ mean(@view(waves[:,end-n_samples:end]), dims = 2)) ./ 2
    for j in axes(waves, 2)
        waves[:,j] -= ms
    end
end
