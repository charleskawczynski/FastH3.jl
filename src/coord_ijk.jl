# IJK hexagon coordinate system functions

const INT32_MAX_3 = typemax(Int32) ÷ Int32(3)

struct CoordIJK
    i::Int32
    j::Int32
    k::Int32
end
CoordIJK() = CoordIJK(Int32(0), Int32(0), Int32(0))
CoordIJK(i::Int, j::Int, k::Int) = CoordIJK(Int32(i), Int32(j), Int32(k))

@enum Direction::Int32 begin
    CENTER_DIGIT = 0
    K_AXES_DIGIT = 1
    J_AXES_DIGIT = 2
    JK_AXES_DIGIT = 3
    I_AXES_DIGIT = 4
    IK_AXES_DIGIT = 5
    IJ_AXES_DIGIT = 6
    INVALID_DIGIT = 7
end

const NUM_DIGITS = INVALID_DIGIT
const PENTAGON_SKIPPED_DIGIT = K_AXES_DIGIT

const UNIT_VECS = [
    CoordIJK(0, 0, 0),  # CENTER_DIGIT (0)
    CoordIJK(0, 0, 1),  # K_AXES_DIGIT (1)
    CoordIJK(0, 1, 0),  # J_AXES_DIGIT (2)
    CoordIJK(0, 1, 1),  # JK_AXES_DIGIT (3)
    CoordIJK(1, 0, 0),  # I_AXES_DIGIT (4)
    CoordIJK(1, 0, 1),  # IK_AXES_DIGIT (5)
    CoordIJK(1, 1, 0),  # IJ_AXES_DIGIT (6)
]

function _hex2dToCoordIJK(v::Vec2d)::CoordIJK
    hk = Int32(0)

    a1 = abs(v.x)
    a2 = abs(v.y)

    x2 = a2 * M_RSIN60
    x1 = a1 + x2 / 2.0

    m1 = trunc(Int, x1)
    m2 = trunc(Int, x2)

    r1 = x1 - m1
    r2 = x2 - m2

    if r1 < 0.5
        if r1 < 1.0 / 3.0
            if r2 < (1.0 + r1) / 2.0
                hi = Int32(m1)
                hj = Int32(m2)
            else
                hi = Int32(m1)
                hj = Int32(m2 + 1)
            end
        else
            if r2 < (1.0 - r1)
                hj = Int32(m2)
            else
                hj = Int32(m2 + 1)
            end

            if (1.0 - r1) <= r2 && r2 < (2.0 * r1)
                hi = Int32(m1 + 1)
            else
                hi = Int32(m1)
            end
        end
    else
        if r1 < 2.0 / 3.0
            if r2 < (1.0 - r1)
                hj = Int32(m2)
            else
                hj = Int32(m2 + 1)
            end

            if (2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1)
                hi = Int32(m1)
            else
                hi = Int32(m1 + 1)
            end
        else
            if r2 < (r1 / 2.0)
                hi = Int32(m1 + 1)
                hj = Int32(m2)
            else
                hi = Int32(m1 + 1)
                hj = Int32(m2 + 1)
            end
        end
    end

    # fold across axes if necessary
    if v.x < 0.0
        if (hj % 2) == 0  # even
            axisi = hj ÷ Int32(2)
            diff = hi - axisi
            hi = Int32(hi - 2 * diff)
        else
            axisi = (hj + Int32(1)) ÷ Int32(2)
            diff = hi - axisi
            hi = Int32(hi - (2 * diff + 1))
        end
    end

    if v.y < 0.0
        hi = hi - (2 * hj + Int32(1)) ÷ Int32(2)
        hj = -hj
    end

    return _ijkNormalize(CoordIJK(hi, hj, hk))
end

function _ijkToHex2d(h::CoordIJK)::Vec2d
    i = h.i - h.k
    j = h.j - h.k
    return Vec2d(Float64(i) - 0.5 * Float64(j), Float64(j) * M_SQRT3_2)
end

function _ijkMatches(c1::CoordIJK, c2::CoordIJK)::Bool
    return c1.i == c2.i && c1.j == c2.j && c1.k == c2.k
end

function _ijkAdd(h1::CoordIJK, h2::CoordIJK)::CoordIJK
    return CoordIJK(h1.i + h2.i, h1.j + h2.j, h1.k + h2.k)
end

function _ijkSub(h1::CoordIJK, h2::CoordIJK)::CoordIJK
    return CoordIJK(h1.i - h2.i, h1.j - h2.j, h1.k - h2.k)
end

function _ijkScale(c::CoordIJK, factor::Int32)::CoordIJK
    return CoordIJK(c.i * factor, c.j * factor, c.k * factor)
end
_ijkScale(c::CoordIJK, factor::Int) = _ijkScale(c, Int32(factor))

function _ijkNormalizeCouldOverflow(ijk::CoordIJK)::Bool
    if ijk.i > ijk.j
        mx, mn = ijk.i, ijk.j
    else
        mx, mn = ijk.j, ijk.i
    end
    if mn < 0
        if ADD_INT32S_OVERFLOWS(mx, mn)
            return true
        end
        if SUB_INT32S_OVERFLOWS(Int32(0), mn)
            return true
        end
        if SUB_INT32S_OVERFLOWS(mx, mn)
            return true
        end
    end
    return false
end

function _ijkNormalize(c::CoordIJK)::CoordIJK
    i, j, k = c.i, c.j, c.k
    if i < 0
        j -= i
        k -= i
        i = Int32(0)
    end
    if j < 0
        i -= j
        k -= j
        j = Int32(0)
    end
    if k < 0
        i -= k
        j -= k
        k = Int32(0)
    end
    mn = min(i, j, k)
    if mn > 0
        i -= mn
        j -= mn
        k -= mn
    end
    return CoordIJK(i, j, k)
end

function _unitIjkToDigit(ijk::CoordIJK)::Direction
    c = _ijkNormalize(ijk)
    for i in Int32(0):Int32(6)
        if _ijkMatches(c, UNIT_VECS[i + 1])
            return Direction(i)
        end
    end
    return INVALID_DIGIT
end

function _upAp7Checked(ijk::CoordIJK)::Tuple{H3Error, CoordIJK}
    i = ijk.i - ijk.k
    j = ijk.j - ijk.k

    if i >= INT32_MAX_3 || j >= INT32_MAX_3 || i < 0 || j < 0
        if ADD_INT32S_OVERFLOWS(i, i)
            return (E_FAILED, ijk)
        end
        i2 = i + i
        if ADD_INT32S_OVERFLOWS(i2, i)
            return (E_FAILED, ijk)
        end
        if ADD_INT32S_OVERFLOWS(j, j)
            return (E_FAILED, ijk)
        end
        j2 = j + j
        if SUB_INT32S_OVERFLOWS(i2 + i, j)
            return (E_FAILED, ijk)
        end
        if ADD_INT32S_OVERFLOWS(i, j2)
            return (E_FAILED, ijk)
        end
    end

    new_i = Int32(round(Int, (i * 3 - j) * M_ONESEVENTH))
    new_j = Int32(round(Int, (i + j * 2) * M_ONESEVENTH))
    result = CoordIJK(new_i, new_j, Int32(0))

    if _ijkNormalizeCouldOverflow(result)
        return (E_FAILED, result)
    end
    return (E_SUCCESS, _ijkNormalize(result))
end

function _upAp7rChecked(ijk::CoordIJK)::Tuple{H3Error, CoordIJK}
    i = ijk.i - ijk.k
    j = ijk.j - ijk.k

    if i >= INT32_MAX_3 || j >= INT32_MAX_3 || i < 0 || j < 0
        if ADD_INT32S_OVERFLOWS(i, i)
            return (E_FAILED, ijk)
        end
        if ADD_INT32S_OVERFLOWS(j, j)
            return (E_FAILED, ijk)
        end
        j2 = j + j
        if ADD_INT32S_OVERFLOWS(j2, j)
            return (E_FAILED, ijk)
        end
        i2 = i + i
        if ADD_INT32S_OVERFLOWS(i2, j)
            return (E_FAILED, ijk)
        end
        if SUB_INT32S_OVERFLOWS(j2 + j, i)
            return (E_FAILED, ijk)
        end
    end

    new_i = Int32(round(Int, (i * 2 + j) * M_ONESEVENTH))
    new_j = Int32(round(Int, (j * 3 - i) * M_ONESEVENTH))
    result = CoordIJK(new_i, new_j, Int32(0))

    if _ijkNormalizeCouldOverflow(result)
        return (E_FAILED, result)
    end
    return (E_SUCCESS, _ijkNormalize(result))
end

function _upAp7(ijk::CoordIJK)::CoordIJK
    i = ijk.i - ijk.k
    j = ijk.j - ijk.k
    new_i = Int32(round(Int, (3 * i - j) * M_ONESEVENTH))
    new_j = Int32(round(Int, (i + 2 * j) * M_ONESEVENTH))
    return _ijkNormalize(CoordIJK(new_i, new_j, Int32(0)))
end

function _upAp7r(ijk::CoordIJK)::CoordIJK
    i = ijk.i - ijk.k
    j = ijk.j - ijk.k
    new_i = Int32(round(Int, (2 * i + j) * M_ONESEVENTH))
    new_j = Int32(round(Int, (3 * j - i) * M_ONESEVENTH))
    return _ijkNormalize(CoordIJK(new_i, new_j, Int32(0)))
end

function _downAp7(ijk::CoordIJK)::CoordIJK
    iVec = _ijkScale(CoordIJK(3, 0, 1), ijk.i)
    jVec = _ijkScale(CoordIJK(1, 3, 0), ijk.j)
    kVec = _ijkScale(CoordIJK(0, 1, 3), ijk.k)
    return _ijkNormalize(_ijkAdd(_ijkAdd(iVec, jVec), kVec))
end

function _downAp7r(ijk::CoordIJK)::CoordIJK
    iVec = _ijkScale(CoordIJK(3, 1, 0), ijk.i)
    jVec = _ijkScale(CoordIJK(0, 3, 1), ijk.j)
    kVec = _ijkScale(CoordIJK(1, 0, 3), ijk.k)
    return _ijkNormalize(_ijkAdd(_ijkAdd(iVec, jVec), kVec))
end

function _neighbor(ijk::CoordIJK, digit::Direction)::CoordIJK
    if digit > CENTER_DIGIT && digit < NUM_DIGITS
        uv = UNIT_VECS[Int(digit) + 1]
        return _ijkNormalize(CoordIJK(ijk.i + uv.i, ijk.j + uv.j, ijk.k + uv.k))
    end
    return ijk
end

function _ijkRotate60ccw(ijk::CoordIJK)::CoordIJK
    iVec = _ijkScale(CoordIJK(1, 1, 0), ijk.i)
    jVec = _ijkScale(CoordIJK(0, 1, 1), ijk.j)
    kVec = _ijkScale(CoordIJK(1, 0, 1), ijk.k)
    return _ijkNormalize(_ijkAdd(_ijkAdd(iVec, jVec), kVec))
end

function _ijkRotate60cw(ijk::CoordIJK)::CoordIJK
    iVec = _ijkScale(CoordIJK(1, 0, 1), ijk.i)
    jVec = _ijkScale(CoordIJK(1, 1, 0), ijk.j)
    kVec = _ijkScale(CoordIJK(0, 1, 1), ijk.k)
    return _ijkNormalize(_ijkAdd(_ijkAdd(iVec, jVec), kVec))
end

function _rotate60ccw(digit::Direction)::Direction
    digit == K_AXES_DIGIT && return IK_AXES_DIGIT
    digit == IK_AXES_DIGIT && return I_AXES_DIGIT
    digit == I_AXES_DIGIT && return IJ_AXES_DIGIT
    digit == IJ_AXES_DIGIT && return J_AXES_DIGIT
    digit == J_AXES_DIGIT && return JK_AXES_DIGIT
    digit == JK_AXES_DIGIT && return K_AXES_DIGIT
    return digit
end

function _rotate60cw(digit::Direction)::Direction
    digit == K_AXES_DIGIT && return JK_AXES_DIGIT
    digit == JK_AXES_DIGIT && return J_AXES_DIGIT
    digit == J_AXES_DIGIT && return IJ_AXES_DIGIT
    digit == IJ_AXES_DIGIT && return I_AXES_DIGIT
    digit == I_AXES_DIGIT && return IK_AXES_DIGIT
    digit == IK_AXES_DIGIT && return K_AXES_DIGIT
    return digit
end

function _downAp3(ijk::CoordIJK)::CoordIJK
    iVec = _ijkScale(CoordIJK(2, 0, 1), ijk.i)
    jVec = _ijkScale(CoordIJK(1, 2, 0), ijk.j)
    kVec = _ijkScale(CoordIJK(0, 1, 2), ijk.k)
    return _ijkNormalize(_ijkAdd(_ijkAdd(iVec, jVec), kVec))
end

function _downAp3r(ijk::CoordIJK)::CoordIJK
    iVec = _ijkScale(CoordIJK(2, 1, 0), ijk.i)
    jVec = _ijkScale(CoordIJK(0, 2, 1), ijk.j)
    kVec = _ijkScale(CoordIJK(1, 0, 2), ijk.k)
    return _ijkNormalize(_ijkAdd(_ijkAdd(iVec, jVec), kVec))
end

function ijkDistance(c1::CoordIJK, c2::CoordIJK)::Int
    diff = _ijkNormalize(_ijkSub(c1, c2))
    return max(abs(diff.i), abs(diff.j), abs(diff.k))
end

function ijkToIj(ijk::CoordIJK)::CoordIJ
    return CoordIJ(Int32(ijk.i - ijk.k), Int32(ijk.j - ijk.k))
end

function ijToIjk(ij::CoordIJ)::Tuple{H3Error, CoordIJK}
    ijk = CoordIJK(ij.i, ij.j, Int32(0))
    if _ijkNormalizeCouldOverflow(ijk)
        return (E_FAILED, ijk)
    end
    return (E_SUCCESS, _ijkNormalize(ijk))
end

function ijkToCube(ijk::CoordIJK)::CoordIJK
    new_i = -ijk.i + ijk.k
    new_j = ijk.j - ijk.k
    new_k = -new_i - new_j
    return CoordIJK(new_i, new_j, new_k)
end

function cubeToIjk(ijk::CoordIJK)::CoordIJK
    return _ijkNormalize(CoordIJK(-ijk.i, ijk.j, Int32(0)))
end
