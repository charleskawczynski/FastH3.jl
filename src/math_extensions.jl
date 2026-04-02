# Math extension functions

"""
Integer exponentiation via binary exponentiation.
"""
function _ipow(base::Int64, exp::Int64)::Int64
    result = Int64(1)
    b = base
    e = exp
    while e != 0
        if (e & 1) != 0
            result *= b
        end
        e >>= 1
        b *= b
    end
    return result
end

function ADD_INT32S_OVERFLOWS(a::Int32, b::Int32)::Bool
    if a > 0
        return typemax(Int32) - a < b
    else
        return typemin(Int32) - a > b
    end
end

function SUB_INT32S_OVERFLOWS(a::Int32, b::Int32)::Bool
    if a >= 0
        return typemin(Int32) + a >= b
    else
        return typemax(Int32) + a + Int32(1) < b
    end
end
