#########################################################################################
#
# Turn a Vector of Vectors into a single Vector containing all elements of the "embedded"
# Vectors
#
##########################################################################################

function unpackvector(vec::Vector{Vector{T}}) where {T}
    outlen::Int64 = 0
    for v in vec
        outlen += length(v)
    end
    out::Vector{Int64} = Vector{Int64}(undef, outlen)

    pos::Int = 1
    for v in vec
        veclen = length(v)
        out[pos:veclen + pos - 1] = v
        pos += veclen
    end
    return out
end
