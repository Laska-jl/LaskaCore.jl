#--------------------
# Opposite of isempty
#--------------------


function isntempty(vec::Vector{T}) where {T}
    return !isempty(vec)
end
