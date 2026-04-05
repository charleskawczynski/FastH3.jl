# Grid disk algorithms (gridDisk, gridRing, etc.)
# Translated from H3 library algos.c

const DIRECTIONS = Direction[
    J_AXES_DIGIT, JK_AXES_DIGIT, K_AXES_DIGIT,
    IK_AXES_DIGIT, I_AXES_DIGIT, IJ_AXES_DIGIT,
]

const NEXT_RING_DIRECTION = I_AXES_DIGIT

# Lookup tables for neighbor traversal.
# Indexed as [oldDigit + 1, direction + 1] (0-based Direction → 1-based Julia).
# Tables transcribed from H3 C source (algos.c).

# Class II: current digit × direction → new digit
const NEW_DIGIT_II = Direction[
    CENTER_DIGIT  K_AXES_DIGIT  J_AXES_DIGIT  JK_AXES_DIGIT I_AXES_DIGIT  IK_AXES_DIGIT IJ_AXES_DIGIT
    K_AXES_DIGIT  I_AXES_DIGIT  JK_AXES_DIGIT IJ_AXES_DIGIT IK_AXES_DIGIT J_AXES_DIGIT  CENTER_DIGIT
    J_AXES_DIGIT  JK_AXES_DIGIT K_AXES_DIGIT  I_AXES_DIGIT  IJ_AXES_DIGIT CENTER_DIGIT  IK_AXES_DIGIT
    JK_AXES_DIGIT IJ_AXES_DIGIT I_AXES_DIGIT  IK_AXES_DIGIT CENTER_DIGIT  K_AXES_DIGIT  J_AXES_DIGIT
    I_AXES_DIGIT  IK_AXES_DIGIT IJ_AXES_DIGIT CENTER_DIGIT  J_AXES_DIGIT  JK_AXES_DIGIT K_AXES_DIGIT
    IK_AXES_DIGIT J_AXES_DIGIT  CENTER_DIGIT  K_AXES_DIGIT  JK_AXES_DIGIT IJ_AXES_DIGIT I_AXES_DIGIT
    IJ_AXES_DIGIT CENTER_DIGIT  IK_AXES_DIGIT J_AXES_DIGIT  K_AXES_DIGIT  I_AXES_DIGIT  JK_AXES_DIGIT
]

# Class II: current digit × direction → coarser-level carry direction
const NEW_ADJUSTMENT_II = Direction[
    CENTER_DIGIT CENTER_DIGIT   CENTER_DIGIT   CENTER_DIGIT    CENTER_DIGIT  CENTER_DIGIT    CENTER_DIGIT
    CENTER_DIGIT K_AXES_DIGIT   CENTER_DIGIT   K_AXES_DIGIT    CENTER_DIGIT  IK_AXES_DIGIT   CENTER_DIGIT
    CENTER_DIGIT CENTER_DIGIT   J_AXES_DIGIT   JK_AXES_DIGIT   CENTER_DIGIT  CENTER_DIGIT    J_AXES_DIGIT
    CENTER_DIGIT K_AXES_DIGIT   JK_AXES_DIGIT  JK_AXES_DIGIT   CENTER_DIGIT  CENTER_DIGIT    CENTER_DIGIT
    CENTER_DIGIT CENTER_DIGIT   CENTER_DIGIT   CENTER_DIGIT    I_AXES_DIGIT  I_AXES_DIGIT    IJ_AXES_DIGIT
    CENTER_DIGIT IK_AXES_DIGIT  CENTER_DIGIT   CENTER_DIGIT    I_AXES_DIGIT  IK_AXES_DIGIT   CENTER_DIGIT
    CENTER_DIGIT CENTER_DIGIT   J_AXES_DIGIT   CENTER_DIGIT    IJ_AXES_DIGIT CENTER_DIGIT    IJ_AXES_DIGIT
]

# Class III: current digit × direction → new digit (cyclic: (i+j) mod 7)
const NEW_DIGIT_III = Direction[
    CENTER_DIGIT  K_AXES_DIGIT  J_AXES_DIGIT  JK_AXES_DIGIT I_AXES_DIGIT  IK_AXES_DIGIT IJ_AXES_DIGIT
    K_AXES_DIGIT  J_AXES_DIGIT  JK_AXES_DIGIT I_AXES_DIGIT  IK_AXES_DIGIT IJ_AXES_DIGIT CENTER_DIGIT
    J_AXES_DIGIT  JK_AXES_DIGIT I_AXES_DIGIT  IK_AXES_DIGIT IJ_AXES_DIGIT CENTER_DIGIT  K_AXES_DIGIT
    JK_AXES_DIGIT I_AXES_DIGIT  IK_AXES_DIGIT IJ_AXES_DIGIT CENTER_DIGIT  K_AXES_DIGIT  J_AXES_DIGIT
    I_AXES_DIGIT  IK_AXES_DIGIT IJ_AXES_DIGIT CENTER_DIGIT  K_AXES_DIGIT  J_AXES_DIGIT  JK_AXES_DIGIT
    IK_AXES_DIGIT IJ_AXES_DIGIT CENTER_DIGIT  K_AXES_DIGIT  J_AXES_DIGIT  JK_AXES_DIGIT I_AXES_DIGIT
    IJ_AXES_DIGIT CENTER_DIGIT  K_AXES_DIGIT  J_AXES_DIGIT  JK_AXES_DIGIT I_AXES_DIGIT  IK_AXES_DIGIT
]

# Class III: current digit × direction → coarser-level carry direction
const NEW_ADJUSTMENT_III = Direction[
    CENTER_DIGIT CENTER_DIGIT    CENTER_DIGIT   CENTER_DIGIT    CENTER_DIGIT   CENTER_DIGIT    CENTER_DIGIT
    CENTER_DIGIT K_AXES_DIGIT    CENTER_DIGIT   JK_AXES_DIGIT   CENTER_DIGIT   K_AXES_DIGIT    CENTER_DIGIT
    CENTER_DIGIT CENTER_DIGIT    J_AXES_DIGIT   J_AXES_DIGIT    CENTER_DIGIT   CENTER_DIGIT    IJ_AXES_DIGIT
    CENTER_DIGIT JK_AXES_DIGIT   J_AXES_DIGIT   JK_AXES_DIGIT   CENTER_DIGIT   CENTER_DIGIT    CENTER_DIGIT
    CENTER_DIGIT CENTER_DIGIT    CENTER_DIGIT   CENTER_DIGIT    I_AXES_DIGIT   IK_AXES_DIGIT   I_AXES_DIGIT
    CENTER_DIGIT K_AXES_DIGIT    CENTER_DIGIT   CENTER_DIGIT    IK_AXES_DIGIT  IK_AXES_DIGIT   CENTER_DIGIT
    CENTER_DIGIT CENTER_DIGIT    IJ_AXES_DIGIT  CENTER_DIGIT    I_AXES_DIGIT   CENTER_DIGIT    IJ_AXES_DIGIT
]

"""
    h3NeighborRotations(origin::H3Index, dir::Direction, rotations::Int) -> (H3Error, H3Index, Int)

Compute the neighbor of `origin` in direction `dir`, applying `rotations` initial
CCW 60° rotations. Returns the error code, the neighbor cell, and the updated
rotation count.

Translated from C source: [`algos.c`](https://github.com/uber/h3/blob/master/src/h3lib/lib/algos.c)
"""
function h3NeighborRotations(origin::H3Index, dir::Direction, rotations::Int)::Tuple{H3Error, H3Index, Int}
    if dir >= INVALID_DIGIT
        return (E_FAILED, H3_NULL, 0)
    end

    current = origin
    d = dir

    rotations = rotations % 6
    for _ in 1:rotations
        d = _rotate60ccw(d)
    end

    newRotations = 0
    oldBaseCell = h3_get_base_cell(current)
    if oldBaseCell < 0 || oldBaseCell >= NUM_BASE_CELLS
        return (E_CELL_INVALID, H3_NULL, 0)
    end
    oldLeadingDigit = _h3LeadingNonZeroDigit(current)

    r = h3_get_resolution(current) - 1
    while true
        if r == -1
            current = h3_set_base_cell(current, _getBaseCellNeighbor(oldBaseCell, d))
            newRotations = Int(baseCellNeighbor60CCWRots[oldBaseCell + 1][Int(d) + 1])

            if h3_get_base_cell(current) == INVALID_BASE_CELL
                # Adjust for the deleted k vertex at the base cell level.
                current = h3_set_base_cell(current, _getBaseCellNeighbor(oldBaseCell, IK_AXES_DIGIT))
                newRotations = Int(baseCellNeighbor60CCWRots[oldBaseCell + 1][Int(IK_AXES_DIGIT) + 1])
                current = _h3Rotate60ccw(current)
                rotations = rotations + 1
            end
            break
        else
            oldDigit = h3_get_index_digit(current, r + 1)
            if oldDigit == INVALID_DIGIT
                return (E_CELL_INVALID, H3_NULL, 0)
            end
            # NOTE: Class III uses _II tables and vice versa (per C source)
            if isResolutionClassIII(r + 1)
                current = h3_set_index_digit(current, r + 1, NEW_DIGIT_II[Int(oldDigit) + 1, Int(d) + 1])
                nextDir = NEW_ADJUSTMENT_II[Int(oldDigit) + 1, Int(d) + 1]
            else
                current = h3_set_index_digit(current, r + 1, NEW_DIGIT_III[Int(oldDigit) + 1, Int(d) + 1])
                nextDir = NEW_ADJUSTMENT_III[Int(oldDigit) + 1, Int(d) + 1]
            end

            if nextDir != CENTER_DIGIT
                d = nextDir
                r -= 1
            else
                break
            end
        end
    end

    newBaseCell = h3_get_base_cell(current)

    if _isBaseCellPentagon(newBaseCell)
        alreadyAdjustedKSubsequence = 0

        if _h3LeadingNonZeroDigit(current) == K_AXES_DIGIT
            if oldBaseCell != newBaseCell
                if _baseCellIsCwOffset(newBaseCell, Int(baseCellData[oldBaseCell + 1].homeFijk.face))
                    current = _h3Rotate60cw(current)
                else
                    current = _h3Rotate60ccw(current)
                end
                alreadyAdjustedKSubsequence = 1
            else
                if oldLeadingDigit == CENTER_DIGIT
                    return (E_PENTAGON, H3_NULL, 0)
                elseif oldLeadingDigit == JK_AXES_DIGIT
                    current = _h3Rotate60ccw(current)
                    rotations = rotations + 1
                elseif oldLeadingDigit == IK_AXES_DIGIT
                    current = _h3Rotate60cw(current)
                    rotations = rotations + 5
                else
                    return (E_FAILED, H3_NULL, 0)
                end
            end
        end

        for _ in 1:newRotations
            current = _h3RotatePent60ccw(current)
        end

        if oldBaseCell != newBaseCell
            if _isBaseCellPolarPentagon(newBaseCell)
                if oldBaseCell != 118 && oldBaseCell != 8 &&
                   _h3LeadingNonZeroDigit(current) != JK_AXES_DIGIT
                    rotations = rotations + 1
                end
            elseif _h3LeadingNonZeroDigit(current) == IK_AXES_DIGIT &&
                   alreadyAdjustedKSubsequence == 0
                rotations = rotations + 1
            end
        end
    else
        for _ in 1:newRotations
            current = _h3Rotate60ccw(current)
        end
    end

    rotations = (rotations + newRotations) % 6
    return (E_SUCCESS, current, rotations)
end

"""
    maxGridDiskSize(k::Int) -> (H3Error, Int64)

Get the maximum number of cells that can be returned by [`gridDisk`](@ref) with radius `k`.
Formula: `3k² + 3k + 1`.

See also the H3 C API: [`maxGridDiskSize`](https://h3geo.org/docs/api/traversal#maxgriddisksize)
"""
function maxGridDiskSize(k::Int)::Tuple{H3Error, Int64}
    if k < 0
        return (E_DOMAIN, Int64(0))
    end
    n = Int64(k)
    return (E_SUCCESS, 3 * n * n + 3 * n + 1)
end

function maxGridRingSize(k::Int)::Tuple{H3Error, Int64}
    if k < 0
        return (E_DOMAIN, Int64(0))
    end
    return (E_SUCCESS, k == 0 ? Int64(1) : Int64(6) * Int64(k))
end

"""
    gridDiskDistancesUnsafe(origin::H3Index, k::Int) -> (H3Error, Vector{H3Index}, Vector{Int})

Produce all cells within grid distance `k` of `origin`, along with their distances.
Fails with `E_PENTAGON` if a pentagon is encountered. Use [`gridDiskDistances`](@ref) for
a safe fallback.

See also the H3 C API: [`gridDiskDistancesUnsafe`](https://h3geo.org/docs/api/traversal#griddiskdistancesunsafe)
"""
function gridDiskDistancesUnsafe(origin::H3Index, k::Int)::Tuple{H3Error, Vector{H3Index}, Vector{Int}}
    err, maxSize = maxGridDiskSize(k)
    if err != E_SUCCESS
        return (err, H3Index[], Int[])
    end

    out = Vector{H3Index}(undef, maxSize)
    distances = Vector{Int}(undef, maxSize)
    idx = 1

    out[idx] = origin
    distances[idx] = 0
    idx += 1

    if isPentagon(origin)
        return (E_PENTAGON, H3Index[], Int[])
    end

    rotations = 0

    for ring in 1:k
        err2, origin, rotations = h3NeighborRotations(origin, NEXT_RING_DIRECTION, rotations)
        if err2 != E_SUCCESS
            return (err2, H3Index[], Int[])
        end

        if isPentagon(origin)
            return (E_PENTAGON, H3Index[], Int[])
        end

        for direction in 1:6
            for _ in 1:ring
                err2, origin, rotations = h3NeighborRotations(origin, DIRECTIONS[direction], rotations)
                if err2 != E_SUCCESS
                    return (err2, H3Index[], Int[])
                end

                out[idx] = origin
                distances[idx] = ring
                idx += 1

                if isPentagon(origin)
                    return (E_PENTAGON, H3Index[], Int[])
                end
            end
        end
    end

    return (E_SUCCESS, out, distances)
end

"""
    gridDiskUnsafe(origin::H3Index, k::Int) -> (H3Error, Vector{H3Index})

Produce all cells within grid distance `k` of `origin`.
Fails with `E_PENTAGON` if a pentagon is encountered. Use [`gridDisk`](@ref) for
a safe fallback.

See also the H3 C API: [`gridDiskUnsafe`](https://h3geo.org/docs/api/traversal#gridDiskunsafe)
"""
function gridDiskUnsafe(origin::H3Index, k::Int)::Tuple{H3Error, Vector{H3Index}}
    err, out, _ = gridDiskDistancesUnsafe(origin, k)
    return (err, out)
end

"""
    gridDiskDistancesSafe(origin::H3Index, k::Int) -> (H3Error, Vector{H3Index}, Vector{Int})

BFS-based safe implementation of [`gridDiskDistances`](@ref) that handles pentagons.
"""
function gridDiskDistancesSafe(origin::H3Index, k::Int)::Tuple{H3Error, Vector{H3Index}, Vector{Int}}
    out = H3Index[]
    distances = Int[]
    visited = Set{H3Index}()

    push!(out, origin)
    push!(distances, 0)
    push!(visited, origin)

    startIdx = 1
    for curDist in 0:(k - 1)
        endIdx = length(out)
        for i in startIdx:endIdx
            cell = out[i]
            for d in 1:6
                err, neighbor, _ = h3NeighborRotations(cell, Direction(d), 0)
                if err == E_PENTAGON
                    continue
                end
                if err != E_SUCCESS
                    return (err, H3Index[], Int[])
                end
                if neighbor ∉ visited
                    push!(visited, neighbor)
                    push!(out, neighbor)
                    push!(distances, curDist + 1)
                end
            end
        end
        startIdx = endIdx + 1
    end

    return (E_SUCCESS, out, distances)
end

"""
    gridDisk(origin::H3Index, k::Int) -> (H3Error, Vector{H3Index})

Produce all cells within grid distance `k` of `origin`. Tries the fast
`gridDiskUnsafe` first and falls back to a safe BFS traversal if pentagons
are encountered.

See also the H3 C API: [`gridDisk`](https://h3geo.org/docs/api/traversal#griddisk)
"""
function gridDisk(origin::H3Index, k::Int)::Tuple{H3Error, Vector{H3Index}}
    err, out = gridDiskUnsafe(origin, k)
    if err == E_SUCCESS
        return (E_SUCCESS, out)
    end
    err2, out2, _ = gridDiskDistancesSafe(origin, k)
    return (err2, out2)
end

"""
    gridDiskDistances(origin::H3Index, k::Int) -> (H3Error, Vector{H3Index}, Vector{Int})

Produce all cells within grid distance `k` of `origin`, along with their distances.
Tries the fast unsafe path first and falls back to safe BFS if needed.

See also the H3 C API: [`gridDiskDistances`](https://h3geo.org/docs/api/traversal#griddiskdistances)
"""
function gridDiskDistances(origin::H3Index, k::Int)::Tuple{H3Error, Vector{H3Index}, Vector{Int}}
    err, out, dist = gridDiskDistancesUnsafe(origin, k)
    if err == E_SUCCESS
        return (E_SUCCESS, out, dist)
    end
    return gridDiskDistancesSafe(origin, k)
end

"""
    gridRingUnsafe(origin::H3Index, k::Int) -> (H3Error, Vector{H3Index})

Produce the cells forming the ring at exactly grid distance `k` from `origin`.
Fails with `E_PENTAGON` if a pentagon is encountered.

See also the H3 C API: [`gridRingUnsafe`](https://h3geo.org/docs/api/traversal#gridringunsafe)
"""
function gridRingUnsafe(origin::H3Index, k::Int)::Tuple{H3Error, Vector{H3Index}}
    if k < 0
        return (E_DOMAIN, H3Index[])
    end
    if k == 0
        return (E_SUCCESS, H3Index[origin])
    end

    if isPentagon(origin)
        return (E_PENTAGON, H3Index[])
    end

    ringSize = 6 * k
    out = Vector{H3Index}(undef, ringSize)

    ringCell = origin
    rotations = 0

    for _ in 1:k
        err, ringCell, rotations = h3NeighborRotations(ringCell, NEXT_RING_DIRECTION, rotations)
        if err != E_SUCCESS
            return (err, H3Index[])
        end
        if isPentagon(ringCell)
            return (E_PENTAGON, H3Index[])
        end
    end

    lastIndex = ringCell
    out[1] = ringCell
    idx = 2

    for direction in 1:6
        for pos in 1:k
            err, ringCell, rotations = h3NeighborRotations(ringCell, DIRECTIONS[direction], rotations)
            if err != E_SUCCESS
                return (err, H3Index[])
            end

            if pos != k || direction != 6
                out[idx] = ringCell
                idx += 1
            end

            if isPentagon(ringCell)
                return (E_PENTAGON, H3Index[])
            end
        end
    end

    if lastIndex != ringCell
        return (E_FAILED, H3Index[])
    end

    return (E_SUCCESS, out)
end

