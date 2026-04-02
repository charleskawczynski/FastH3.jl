# Local IJ coordinate system functions
# Translated from H3 library localij.c

# Pentagon rotation tables: PENTAGON_ROTATIONS[leadingDigit+1, dir+1]
# -1 means invalid. Value is the number of 60° CCW rotations needed.
const PENTAGON_ROTATIONS = Int32[
     0 -1  0  0  0  0  0;
    -1 -1 -1 -1 -1 -1 -1;
     0 -1  0  0  0  1  0;
     0 -1  0  0  1  1  0;
     0 -1  0  5  0  0  0;
     0 -1  5  5  0  0  0;
     0 -1  0  0  0  0  0
]

const PENTAGON_ROTATIONS_REVERSE = Int32[
     0 -1  0  0  0  0  0;
    -1 -1 -1 -1 -1 -1 -1;
     0 -1  0  0  0  5  0;
     0 -1  0  0  5  5  0;
     0 -1  0  1  0  0  0;
     0 -1  1  1  0  0  0;
     0 -1  0  0  0  0  0
]

"""
Find the neighbor of an H3 cell in the given direction, applying rotations
as needed for base cell transitions and pentagon handling.
Returns (error, neighbor_index, new_rotations).
"""
function _h3NeighborRotations(origin::H3Index, dir::Direction, rotations::Int)::Tuple{H3Error, H3Index, Int}
    current = origin

    d = dir
    for _ in 1:rotations
        d = _rotate60ccw(d)
    end

    if d <= CENTER_DIGIT || d >= INVALID_DIGIT
        return (E_FAILED, H3_NULL, 0)
    end

    res = h3_get_resolution(current)
    oldBaseCell = h3_get_base_cell(current)
    oldLeadingDigit = _h3LeadingNonZeroDigit(current)
    newRotations = 0

    r = res - 1
    while true
        if r == -1
            newBC = _getBaseCellNeighbor(oldBaseCell, d)
            newRotations = Int(baseCellNeighbor60CCWRots[oldBaseCell + 1][Int(d) + 1])

            if newBC == INVALID_BASE_CELL
                newBC = _getBaseCellNeighbor(oldBaseCell, _rotate60cw(d))
                newRotations = Int(baseCellNeighbor60CCWRots[oldBaseCell + 1][Int(_rotate60cw(d)) + 1])
                current = _h3Rotate60ccw(current)
                newRotations += 1
            end

            current = h3_set_base_cell(current, newBC)
            break
        else
            digitPos = r + 1
            oldDigit = h3_get_index_digit(current, digitPos)

            sumIjk = _ijkNormalize(_ijkAdd(UNIT_VECS[Int(oldDigit) + 1], UNIT_VECS[Int(d) + 1]))

            carry = sumIjk
            if isResolutionClassIII(digitPos)
                carry = _upAp7(carry)
                center = _downAp7(carry)
            else
                carry = _upAp7r(carry)
                center = _downAp7r(carry)
            end

            diff = _ijkNormalize(_ijkSub(sumIjk, center))

            newDigit = _unitIjkToDigit(diff)
            nextDir = _unitIjkToDigit(carry)

            if newDigit == INVALID_DIGIT
                return (E_FAILED, H3_NULL, 0)
            end

            current = h3_set_index_digit(current, digitPos, newDigit)

            if nextDir == CENTER_DIGIT
                break
            end

            d = nextDir
            r -= 1
        end
    end

    newBC = h3_get_base_cell(current)
    if _isBaseCellPentagon(newBC)
        alreadyAdjustedKSubsequence = false

        if _h3LeadingNonZeroDigit(current) == K_AXES_DIGIT
            if oldBaseCell != newBC
                current = _h3Rotate60ccw(current)
            else
                current = _h3Rotate60cw(current)
            end
            alreadyAdjustedKSubsequence = true
        end

        for _ in 1:newRotations
            current = _h3RotatePent60ccw(current)
        end

        if alreadyAdjustedKSubsequence && _h3LeadingNonZeroDigit(current) == K_AXES_DIGIT
            if oldBaseCell == newBC
                current = _h3Rotate60cw(current)
            else
                current = _h3Rotate60ccw(current)
            end
        end
    else
        for _ in 1:newRotations
            current = _h3Rotate60ccw(current)
        end
    end

    return (E_SUCCESS, current, newRotations)
end

"""
Compute the local IJK coordinates of `h` relative to `origin`.
Both must be at the same resolution.
"""
function _cellToLocalIjk(origin::H3Index, h::H3Index)::Tuple{H3Error, CoordIJK}
    res = h3_get_resolution(origin)
    if res != h3_get_resolution(h)
        return (E_RES_MISMATCH, CoordIJK())
    end

    originBC = h3_get_base_cell(origin)
    hBC = h3_get_base_cell(h)

    if originBC < 0 || originBC >= NUM_BASE_CELLS
        return (E_CELL_INVALID, CoordIJK())
    end
    if hBC < 0 || hBC >= NUM_BASE_CELLS
        return (E_CELL_INVALID, CoordIJK())
    end

    dir = CENTER_DIGIT
    if originBC != hBC
        dir = _getBaseCellDirection(originBC, hBC)
        if dir == INVALID_DIGIT
            return (E_FAILED, CoordIJK())
        end
    end

    originOnPent = _isBaseCellPentagon(originBC)
    indexOnPent = _isBaseCellPentagon(hBC)

    baseCellRotations = 0
    pentagonRotations = 0

    if dir != CENTER_DIGIT
        baseCellRotations = Int(baseCellNeighbor60CCWRots[originBC + 1][Int(dir) + 1])
        if baseCellRotations < 0
            return (E_FAILED, CoordIJK())
        end
    end

    hIjk = CoordIJK(0, 0, 0)
    if dir != CENTER_DIGIT
        hIjk = UNIT_VECS[Int(dir) + 1]

        if originOnPent
            originLeading = _h3LeadingNonZeroDigit(origin)
            pentRot = PENTAGON_ROTATIONS[Int(originLeading) + 1, Int(dir) + 1]
            if pentRot < 0
                return (E_FAILED, CoordIJK())
            end
            for _ in 1:pentRot
                hIjk = _ijkRotate60ccw(hIjk)
            end
        else
            for _ in 1:baseCellRotations
                hIjk = _ijkRotate60ccw(hIjk)
            end
        end

        if indexOnPent
            pentSearchDir = dir
            if originOnPent
                originLeading = _h3LeadingNonZeroDigit(origin)
                pentRot = PENTAGON_ROTATIONS[Int(originLeading) + 1, Int(dir) + 1]
                if pentRot >= 0
                    for _ in 1:pentRot
                        pentSearchDir = _rotate60ccw(pentSearchDir)
                    end
                end
            else
                for _ in 1:baseCellRotations
                    pentSearchDir = _rotate60ccw(pentSearchDir)
                end
            end

            revDir = _getBaseCellDirection(hBC, originBC)
            if revDir == INVALID_DIGIT
                return (E_FAILED, CoordIJK())
            end

            pentagonRotations = PENTAGON_ROTATIONS_REVERSE[Int(revDir) + 1, Int(pentSearchDir) + 1]
            if pentagonRotations < 0
                return (E_FAILED, CoordIJK())
            end
        end
    end

    for r in 1:res
        if isResolutionClassIII(r)
            hIjk = _downAp7(hIjk)
        else
            hIjk = _downAp7r(hIjk)
        end

        digit = h3_get_index_digit(h, r)

        if indexOnPent
            if pentagonRotations > 0 && digit != CENTER_DIGIT
                for _ in 1:pentagonRotations
                    digit = _rotate60ccw(digit)
                end
                if digit == K_AXES_DIGIT
                    pentagonRotations = 0
                    continue
                end
                pentagonRotations = 0
            end
        end

        hIjk = _neighbor(hIjk, digit)
    end

    originIjk = CoordIJK(0, 0, 0)
    for r in 1:res
        if isResolutionClassIII(r)
            originIjk = _downAp7(originIjk)
        else
            originIjk = _downAp7r(originIjk)
        end
        originIjk = _neighbor(originIjk, h3_get_index_digit(origin, r))
    end

    diff = _ijkNormalize(_ijkSub(hIjk, originIjk))
    return (E_SUCCESS, diff)
end

"""
Convert a cell to local IJ coordinates relative to an origin cell.
Both cells must be at the same resolution.
"""
function cellToLocalIj(origin::H3Index, h::H3Index, mode::UInt32)::Tuple{H3Error, CoordIJ}
    err, ijk = _cellToLocalIjk(origin, h)
    if err != E_SUCCESS
        return (err, CoordIJ())
    end
    return (E_SUCCESS, ijkToIj(ijk))
end

"""
Convert local IJ coordinates back to a cell H3Index, using `origin` as the
local coordinate system anchor.
"""
function localIjToCell(origin::H3Index, ij::CoordIJ, mode::UInt32)::Tuple{H3Error, H3Index}
    err, offsetIjk = ijToIjk(ij)
    if err != E_SUCCESS
        return (err, H3_NULL)
    end

    res = h3_get_resolution(origin)
    originBC = h3_get_base_cell(origin)

    originIjk = CoordIJK(0, 0, 0)
    for r in 1:res
        if isResolutionClassIII(r)
            originIjk = _downAp7(originIjk)
        else
            originIjk = _downAp7r(originIjk)
        end
        originIjk = _neighbor(originIjk, h3_get_index_digit(origin, r))
    end

    targetIjk = _ijkNormalize(_ijkAdd(originIjk, offsetIjk))

    h = H3_INIT
    h = h3_set_mode(h, H3_CELL_MODE)
    h = h3_set_resolution(h, res)

    ijk = targetIjk
    for r in res:-1:1
        lastIJK = ijk
        if isResolutionClassIII(r)
            ijk = _upAp7(ijk)
            lastCenter = _downAp7(ijk)
        else
            ijk = _upAp7r(ijk)
            lastCenter = _downAp7r(ijk)
        end
        diff = _ijkNormalize(_ijkSub(lastIJK, lastCenter))
        h = h3_set_index_digit(h, r, _unitIjkToDigit(diff))
    end

    ijk = _ijkNormalize(ijk)
    baseCellDir = _unitIjkToDigit(ijk)

    if baseCellDir == CENTER_DIGIT
        h = h3_set_base_cell(h, originBC)
    elseif baseCellDir == INVALID_DIGIT
        return (E_FAILED, H3_NULL)
    else
        newBC = _getBaseCellNeighbor(originBC, baseCellDir)
        if newBC == INVALID_BASE_CELL
            return (E_FAILED, H3_NULL)
        end
        h = h3_set_base_cell(h, newBC)

        numRots = Int(baseCellNeighbor60CCWRots[originBC + 1][Int(baseCellDir) + 1])
        if numRots < 0
            return (E_FAILED, H3_NULL)
        end

        if _isBaseCellPentagon(newBC)
            if _h3LeadingNonZeroDigit(h) == K_AXES_DIGIT
                if originBC != newBC
                    h = _h3Rotate60ccw(h)
                else
                    h = _h3Rotate60cw(h)
                end
            end
            for _ in 1:numRots
                h = _h3RotatePent60ccw(h)
            end
        else
            for _ in 1:numRots
                h = _h3Rotate60ccw(h)
            end
        end
    end

    return (E_SUCCESS, h)
end

"""
Compute the grid distance between two H3 cells at the same resolution.
"""
function gridDistance(origin::H3Index, h3::H3Index)::Tuple{H3Error, Int64}
    err, ijk = _cellToLocalIjk(origin, h3)
    if err != E_SUCCESS
        err2, ijk2 = _cellToLocalIjk(h3, origin)
        if err2 != E_SUCCESS
            return (err, Int64(0))
        end
        ijk = ijk2
    end
    dist = ijkDistance(ijk, CoordIJK(0, 0, 0))
    return (E_SUCCESS, Int64(dist))
end

"""
Number of cells in a grid path between `start` and `end_`.
"""
function gridPathCellsSize(start::H3Index, end_::H3Index)::Tuple{H3Error, Int64}
    err, dist = gridDistance(start, end_)
    if err != E_SUCCESS
        return (err, Int64(0))
    end
    return (E_SUCCESS, dist + Int64(1))
end

function _cubeRound(fi::Float64, fj::Float64, fk::Float64)::Tuple{Int32, Int32, Int32}
    ri = round(Int32, fi)
    rj = round(Int32, fj)
    rk = round(Int32, fk)

    i_diff = abs(Float64(ri) - fi)
    j_diff = abs(Float64(rj) - fj)
    k_diff = abs(Float64(rk) - fk)

    if i_diff > j_diff && i_diff > k_diff
        ri = Int32(-rj - rk)
    elseif j_diff > k_diff
        rj = Int32(-ri - rk)
    else
        rk = Int32(-ri - rj)
    end

    return (ri, rj, rk)
end

"""
Compute the grid path (line of cells) between `start_` and `end_`.
Uses cube coordinate interpolation for correct hex grid line drawing.
"""
function gridPathCells(start_::H3Index, end_::H3Index)::Tuple{H3Error, Vector{H3Index}}
    err, dist = gridDistance(start_, end_)
    if err != E_SUCCESS
        return (err, H3Index[])
    end

    if dist == 0
        return (E_SUCCESS, [start_])
    end

    err2, endIjk = _cellToLocalIjk(start_, end_)
    if err2 != E_SUCCESS
        return (err2, H3Index[])
    end

    startCube = ijkToCube(CoordIJK(0, 0, 0))
    endCubeCopy = ijkToCube(endIjk)

    iStep = Float64(endCubeCopy.i - startCube.i) / Float64(dist)
    jStep = Float64(endCubeCopy.j - startCube.j) / Float64(dist)
    kStep = Float64(endCubeCopy.k - startCube.k) / Float64(dist)

    numCells = dist + 1
    path = Vector{H3Index}(undef, numCells)

    for n in Int64(0):(numCells - Int64(1))
        ci = Float64(startCube.i) + iStep * Float64(n)
        cj = Float64(startCube.j) + jStep * Float64(n)
        ck = Float64(startCube.k) + kStep * Float64(n)

        ri, rj, rk = _cubeRound(ci, cj, ck)

        ijk = cubeToIjk(CoordIJK(ri, rj, rk))
        ij = ijkToIj(ijk)

        err3, cell = localIjToCell(start_, ij, UInt32(0))
        if err3 != E_SUCCESS
            return (err3, H3Index[])
        end
        path[n + 1] = cell
    end

    return (E_SUCCESS, path)
end
