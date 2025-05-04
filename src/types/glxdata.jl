
"""
    struct NPTraces{T,U<:Real} <: AbstractMatrix{T}


"""
struct NPTraces{T,U<:Real} <: AbstractMatrix{T}
    traces::Matrix{T}
    samplerate::U
    t_start::U
    t_end::U
end

# Interface implementation
Base.size(A::NPTraces) = size(A.traces)
Base.getindex(A::NPTraces, i::Int) = getindex(A.traces, i)

"""
    traces(data::NPTraces)

TBW
"""
traces(data::NPTraces) = data.traces

"""
    time(data::NPTraces)

TBW
"""
time(data::NPTraces) = (data.t_start, data.t_end)

"""
    t_start(data::NPTraces)

TBW
"""
t_start(data::NPTraces) = data.t_start

"""
    t_end(data::NPTraces)

TBW
"""
t_end(data::NPTraces) = data.t_end

