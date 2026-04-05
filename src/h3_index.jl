# H3Index utility functions - bit manipulation and core API

# Bit layout constants
const H3_NUM_BITS = 64
const H3_MAX_OFFSET = 63
const H3_MODE_OFFSET = 59
const H3_BC_OFFSET = 45
const H3_RES_OFFSET = 52
const H3_RESERVED_OFFSET = 56
const H3_PER_DIGIT_OFFSET = 3

const H3_HIGH_BIT_MASK = UInt64(1) << H3_MAX_OFFSET
const H3_MODE_MASK = UInt64(15) << H3_MODE_OFFSET
const H3_BC_MASK = UInt64(127) << H3_BC_OFFSET
const H3_RES_MASK = UInt64(15) << H3_RES_OFFSET
const H3_RESERVED_MASK = UInt64(7) << H3_RESERVED_OFFSET
const H3_DIGIT_MASK = UInt64(7)
const H3_RESERVED_MASK_NEGATIVE = ~H3_RESERVED_MASK

# H3_INIT: mode 0, res 0, base cell 0, all index digits set to 7
const H3_INIT = UInt64(35184372088831)

# Bit manipulation functions (all return new H3Index, no mutation)
h3_get_high_bit(h::H3Index) = Int((h & H3_HIGH_BIT_MASK) >> H3_MAX_OFFSET)
h3_get_mode(h::H3Index) = Int((h & H3_MODE_MASK) >> H3_MODE_OFFSET)
h3_set_mode(h::H3Index, v::Int) = (h & ~H3_MODE_MASK) | (UInt64(v) << H3_MODE_OFFSET)
h3_get_base_cell(h::H3Index) = Int((h & H3_BC_MASK) >> H3_BC_OFFSET)
h3_set_base_cell(h::H3Index, bc::Int) = (h & ~H3_BC_MASK) | (UInt64(bc) << H3_BC_OFFSET)
h3_get_resolution(h::H3Index) = Int((h & H3_RES_MASK) >> H3_RES_OFFSET)
h3_set_resolution(h::H3Index, res::Int) = (h & ~H3_RES_MASK) | (UInt64(res) << H3_RES_OFFSET)

function h3_get_index_digit(h::H3Index, res::Int)::Direction
    shift = (MAX_H3_RES - res) * H3_PER_DIGIT_OFFSET
    return Direction(Int((h >> shift) & H3_DIGIT_MASK))
end

function h3_set_index_digit(h::H3Index, res::Int, digit)::H3Index
    shift = (MAX_H3_RES - res) * H3_PER_DIGIT_OFFSET
    return (h & ~(H3_DIGIT_MASK << shift)) | (UInt64(Int(digit)) << shift)
end

h3_get_reserved_bits(h::H3Index) = Int((h & H3_RESERVED_MASK) >> H3_RESERVED_OFFSET)
h3_set_reserved_bits(h::H3Index, v::Int) = (h & ~H3_RESERVED_MASK) | (UInt64(v) << H3_RESERVED_OFFSET)

function setH3Index(res::Int, baseCell::Int, initDigit::Direction)::H3Index
    h = H3_INIT
    h = h3_set_mode(h, H3_CELL_MODE)
    h = h3_set_resolution(h, res)
    h = h3_set_base_cell(h, baseCell)
    for r in 1:res
        h = h3_set_index_digit(h, r, initDigit)
    end
    return h
end
setH3Index(res::Int, baseCell::Int, initDigit::Int) = setH3Index(res, baseCell, Direction(initDigit))

"""
    getResolution(h::H3Index) -> Int

Get the resolution of an H3 index.

See also the H3 C API: [`getResolution`](https://h3geo.org/docs/api/inspection#getresolution)
"""
getResolution(h::H3Index)::Int = h3_get_resolution(h)

"""
    getBaseCellNumber(h::H3Index) -> Int

Get the base cell number (0–121) of an H3 index.

See also the H3 C API: [`getBaseCellNumber`](https://h3geo.org/docs/api/inspection#getbasecellnumber)
"""
getBaseCellNumber(h::H3Index)::Int = h3_get_base_cell(h)

"""
    getIndexDigit(h::H3Index, res::Int) -> Tuple{H3Error, Int}

Get the index digit at resolution `res` (1-indexed). Returns `(error, digit)`.
"""
function getIndexDigit(h::H3Index, res::Int)::Tuple{H3Error, Int}
    if res < 1 || res > MAX_H3_RES
        return (E_RES_DOMAIN, 0)
    end
    return (E_SUCCESS, Int(h3_get_index_digit(h, res)))
end

function _h3LeadingNonZeroDigit(h::H3Index)::Direction
    for r in 1:h3_get_resolution(h)
        d = h3_get_index_digit(h, r)
        if d != CENTER_DIGIT
            return d
        end
    end
    return CENTER_DIGIT
end

# Validation helpers
function _hasGoodTopBits(h::H3Index)::Bool
    return (h >> (64 - 8)) == 0b00001000
end

function _hasAny7UptoRes(h::H3Index, res::Int)::Bool
    MHI = UInt64(0b100100100100100100100100100100100100100100100)
    MLO = MHI >> 2
    shift = 3 * (15 - res)
    h >>= shift
    h <<= shift
    h = (h & MHI & (~h - MLO))
    return h != 0
end

function _hasAll7AfterRes(h::H3Index, res::Int)::Bool
    if res < 15
        shift = 19 + 3 * res
        h = ~h
        h <<= shift
        h >>= shift
        return h == 0
    end
    return true
end

function _hasDeletedSubsequence(h::H3Index, base_cell::Int)::Bool
    if _isBaseCellPentagon(base_cell)
        h <<= 19
        h >>= 19
        if h == 0
            return false
        end
        return (63 - leading_zeros(h)) % 3 == 0
    end
    return false
end

"""
    isValidCell(h::H3Index) -> Bool

Check whether an H3 index represents a valid cell.

See also the H3 C API: [`isValidCell`](https://h3geo.org/docs/api/inspection#isvalidcell)
"""
function isValidCell(h::H3Index)::Bool
    _hasGoodTopBits(h) || return false
    res = h3_get_resolution(h)
    bc = h3_get_base_cell(h)
    bc >= NUM_BASE_CELLS && return false
    _hasAny7UptoRes(h, res) && return false
    _hasAll7AfterRes(h, res) || return false
    _hasDeletedSubsequence(h, bc) && return false
    return true
end

"""
    isPentagon(h::H3Index) -> Bool

Check whether an H3 index represents a pentagon cell.

See also the H3 C API: [`isPentagon`](https://h3geo.org/docs/api/inspection#ispentagon)
"""
function isPentagon(h::H3Index)::Bool
    return _isBaseCellPentagon(h3_get_base_cell(h)) &&
           _h3LeadingNonZeroDigit(h) == CENTER_DIGIT
end

"""
    isResClassIII(h::H3Index) -> Bool

Check whether the resolution of an H3 index is Class III.

See also the H3 C API: [`isResClassIII`](https://h3geo.org/docs/api/inspection#isresclassiii)
"""
isResClassIII(h::H3Index)::Bool = h3_get_resolution(h) % 2 != 0

# Rotation functions
function _h3Rotate60ccw(h::H3Index)::H3Index
    res = h3_get_resolution(h)
    for r in 1:res
        oldDigit = h3_get_index_digit(h, r)
        h = h3_set_index_digit(h, r, _rotate60ccw(oldDigit))
    end
    return h
end

function _h3Rotate60cw(h::H3Index)::H3Index
    res = h3_get_resolution(h)
    for r in 1:res
        h = h3_set_index_digit(h, r, _rotate60cw(h3_get_index_digit(h, r)))
    end
    return h
end

function _h3RotatePent60ccw(h::H3Index)::H3Index
    foundFirstNonZeroDigit = false
    res = h3_get_resolution(h)
    for r in 1:res
        h = h3_set_index_digit(h, r, _rotate60ccw(h3_get_index_digit(h, r)))
        if !foundFirstNonZeroDigit && h3_get_index_digit(h, r) != CENTER_DIGIT
            foundFirstNonZeroDigit = true
            if _h3LeadingNonZeroDigit(h) == K_AXES_DIGIT
                h = _h3Rotate60ccw(h)
            end
        end
    end
    return h
end

function _h3RotatePent60cw(h::H3Index)::H3Index
    foundFirstNonZeroDigit = false
    res = h3_get_resolution(h)
    for r in 1:res
        h = h3_set_index_digit(h, r, _rotate60cw(h3_get_index_digit(h, r)))
        if !foundFirstNonZeroDigit && h3_get_index_digit(h, r) != CENTER_DIGIT
            foundFirstNonZeroDigit = true
            if _h3LeadingNonZeroDigit(h) == K_AXES_DIGIT
                h = _h3Rotate60cw(h)
            end
        end
    end
    return h
end

# FaceIJK ↔ H3Index conversions
function _faceIjkToH3(fijk::FaceIJK, res::Int)::H3Index
    h = H3_INIT
    h = h3_set_mode(h, H3_CELL_MODE)
    h = h3_set_resolution(h, res)

    if res == 0
        if fijk.coord.i > MAX_FACE_COORD || fijk.coord.j > MAX_FACE_COORD ||
           fijk.coord.k > MAX_FACE_COORD
            return H3_NULL
        end
        h = h3_set_base_cell(h, _faceIjkToBaseCell(fijk))
        return h
    end

    ijk = fijk.coord

    for r in (res - 1):-1:0
        lastIJK = ijk
        if isResolutionClassIII(r + 1)
            ijk = _upAp7(ijk)
            lastCenter = _downAp7(ijk)
        else
            ijk = _upAp7r(ijk)
            lastCenter = _downAp7r(ijk)
        end
        diff = _ijkNormalize(_ijkSub(lastIJK, lastCenter))
        h = h3_set_index_digit(h, r + 1, _unitIjkToDigit(diff))
    end

    fijkBC = FaceIJK(fijk.face, ijk)

    if ijk.i > MAX_FACE_COORD || ijk.j > MAX_FACE_COORD ||
       ijk.k > MAX_FACE_COORD
        return H3_NULL
    end

    baseCell = _faceIjkToBaseCell(fijkBC)
    h = h3_set_base_cell(h, baseCell)

    numRots = _faceIjkToBaseCellCCWrot60(fijkBC)
    if _isBaseCellPentagon(baseCell)
        if _h3LeadingNonZeroDigit(h) == K_AXES_DIGIT
            if _baseCellIsCwOffset(baseCell, Int(fijkBC.face))
                h = _h3Rotate60cw(h)
            else
                h = _h3Rotate60ccw(h)
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

    return h
end

function _h3ToFaceIjkWithInitializedFijk(h::H3Index, fijk::FaceIJK)::Tuple{Int, FaceIJK}
    ijk = fijk.coord
    res = h3_get_resolution(h)

    possibleOverage = 1
    if !_isBaseCellPentagon(h3_get_base_cell(h)) &&
       (res == 0 || (fijk.coord.i == 0 && fijk.coord.j == 0 && fijk.coord.k == 0))
        possibleOverage = 0
    end

    for r in 1:res
        if isResolutionClassIII(r)
            ijk = _downAp7(ijk)
        else
            ijk = _downAp7r(ijk)
        end
        ijk = _neighbor(ijk, h3_get_index_digit(h, r))
    end

    return (possibleOverage, FaceIJK(fijk.face, ijk))
end

function _h3ToFaceIjk(h::H3Index)::Tuple{H3Error, FaceIJK}
    baseCell = h3_get_base_cell(h)
    if baseCell < 0 || baseCell >= NUM_BASE_CELLS
        return (E_CELL_INVALID, FaceIJK())
    end

    if _isBaseCellPentagon(baseCell) && _h3LeadingNonZeroDigit(h) == IK_AXES_DIGIT
        h = _h3Rotate60cw(h)
    end

    fijk = _baseCellToFaceIjk(baseCell)
    possibleOverage, fijk = _h3ToFaceIjkWithInitializedFijk(h, fijk)
    if possibleOverage == 0
        return (E_SUCCESS, fijk)
    end

    origIJK = fijk.coord
    res = h3_get_resolution(h)

    coord = fijk.coord
    if isResolutionClassIII(res)
        coord = _downAp7r(coord)
        res += 1
    end
    fijk = FaceIJK(fijk.face, coord)

    pentLeading4 = _isBaseCellPentagon(baseCell) && _h3LeadingNonZeroDigit(h) == I_AXES_DIGIT ? 1 : 0
    overage, fijk = _adjustOverageClassII(fijk, res, pentLeading4, 0)
    if overage != NO_OVERAGE
        if _isBaseCellPentagon(baseCell)
            while true
                overage, fijk = _adjustOverageClassII(fijk, res, 0, 0)
                overage == NO_OVERAGE && break
            end
        end
        if res != h3_get_resolution(h)
            fijk = FaceIJK(fijk.face, _upAp7r(fijk.coord))
        end
    elseif res != h3_get_resolution(h)
        fijk = FaceIJK(fijk.face, origIJK)
    end
    return (E_SUCCESS, fijk)
end

# Public API functions
"""
    latLngToCell(g::LatLng, res::Int) -> Tuple{H3Error, H3Index}

Convert a latitude/longitude pair to the containing H3 cell at the given resolution.

See also the H3 C API: [`latLngToCell`](https://h3geo.org/docs/api/indexing#latlngtocell)
"""
function latLngToCell(g::LatLng, res::Int)::Tuple{H3Error, H3Index}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, H3_NULL)
    end
    if !isfinite(g.lat) || !isfinite(g.lng)
        return (E_LATLNG_DOMAIN, H3_NULL)
    end
    fijk = _geoToFaceIjk(g, res)
    out = _faceIjkToH3(fijk, res)
    if out != H3_NULL
        return (E_SUCCESS, out)
    else
        return (E_FAILED, H3_NULL)
    end
end

"""
    cellToLatLng(h::H3Index) -> Tuple{H3Error, LatLng}

Convert an H3 cell to its center latitude/longitude.

See also the H3 C API: [`cellToLatLng`](https://h3geo.org/docs/api/indexing#celltolatlng)
"""
function cellToLatLng(h::H3Index)::Tuple{H3Error, LatLng}
    err, fijk = _h3ToFaceIjk(h)
    if err != E_SUCCESS
        return (err, LatLng())
    end
    g = _faceIjkToGeo(fijk, h3_get_resolution(h))
    return (E_SUCCESS, g)
end

"""
    cellToBoundary(h::H3Index) -> Tuple{H3Error, CellBoundary}

Compute the cell boundary (vertices) for an H3 cell.

See also the H3 C API: [`cellToBoundary`](https://h3geo.org/docs/api/indexing#celltoboundary)
"""
function cellToBoundary(h::H3Index)::Tuple{H3Error, CellBoundary}
    err, fijk = _h3ToFaceIjk(h)
    if err != E_SUCCESS
        return (err, CellBoundary())
    end
    if isPentagon(h)
        cb = _faceIjkPentToCellBoundary(fijk, h3_get_resolution(h), 0, NUM_PENT_VERTS)
    else
        cb = _faceIjkToCellBoundary(fijk, h3_get_resolution(h), 0, NUM_HEX_VERTS)
    end
    return (E_SUCCESS, cb)
end

"""
    cellToParent(h::H3Index, parentRes::Int) -> Tuple{H3Error, H3Index}

Get the parent cell at the given coarser resolution.

See also the H3 C API: [`cellToParent`](https://h3geo.org/docs/api/hierarchy#celltoparent)
"""
function cellToParent(h::H3Index, parentRes::Int)::Tuple{H3Error, H3Index}
    childRes = h3_get_resolution(h)
    if parentRes < 0 || parentRes > MAX_H3_RES
        return (E_RES_DOMAIN, H3_NULL)
    end
    if parentRes > childRes
        return (E_RES_MISMATCH, H3_NULL)
    end
    if parentRes == childRes
        return (E_SUCCESS, h)
    end
    parentH = h3_set_resolution(h, parentRes)
    for i in (parentRes + 1):childRes
        parentH = h3_set_index_digit(parentH, i, H3_DIGIT_MASK)
    end
    return (E_SUCCESS, parentH)
end

"""
    cellToChildrenSize(h::H3Index, childRes::Int) -> Tuple{H3Error, Int64}

Get the number of children at the given finer resolution.

See also the H3 C API: [`cellToChildrenSize`](https://h3geo.org/docs/api/hierarchy#celltochildrensize)
"""
function cellToChildrenSize(h::H3Index, childRes::Int)::Tuple{H3Error, Int64}
    parentRes = h3_get_resolution(h)
    if childRes < parentRes || childRes > MAX_H3_RES
        return (E_RES_DOMAIN, Int64(0))
    end
    n = childRes - parentRes
    if isPentagon(h)
        return (E_SUCCESS, Int64(1) + Int64(5) * (_ipow(Int64(7), Int64(n)) - Int64(1)) ÷ Int64(6))
    else
        return (E_SUCCESS, _ipow(Int64(7), Int64(n)))
    end
end

function makeDirectChild(h::H3Index, cellNumber::Int)::H3Index
    childRes = h3_get_resolution(h) + 1
    childH = h3_set_resolution(h, childRes)
    childH = h3_set_index_digit(childH, childRes, cellNumber)
    return childH
end

function _zeroIndexDigits(h::H3Index, start::Int, end_::Int)::H3Index
    if start > end_
        return h
    end
    m = ~UInt64(0)
    m <<= H3_PER_DIGIT_OFFSET * (end_ - start + 1)
    m = ~m
    m <<= H3_PER_DIGIT_OFFSET * (MAX_H3_RES - end_)
    m = ~m
    return h & m
end

"""
    cellToCenterChild(h::H3Index, childRes::Int) -> Tuple{H3Error, H3Index}

Get the center child cell at the given finer resolution.

See also the H3 C API: [`cellToCenterChild`](https://h3geo.org/docs/api/hierarchy#celltocenterchild)
"""
function cellToCenterChild(h::H3Index, childRes::Int)::Tuple{H3Error, H3Index}
    parentRes = h3_get_resolution(h)
    if childRes < parentRes || childRes > MAX_H3_RES
        return (E_RES_DOMAIN, H3_NULL)
    end
    h = _zeroIndexDigits(h, parentRes + 1, childRes)
    h = h3_set_resolution(h, childRes)
    return (E_SUCCESS, h)
end

"""
    maxFaceCount(h::H3Index) -> Tuple{H3Error, Int}

Get the maximum number of icosahedron faces a cell may intersect.

See also the H3 C API: [`maxFaceCount`](https://h3geo.org/docs/api/inspection#maxfacecount)
"""
function maxFaceCount(h::H3Index)::Tuple{H3Error, Int}
    return (E_SUCCESS, isPentagon(h) ? 5 : 2)
end

"""
    pentagonCount() -> Int

Return the number of pentagon cells per resolution (always 12).

See also the H3 C API: [`pentagonCount`](https://h3geo.org/docs/api/misc#pentagoncount)
"""
function pentagonCount()::Int
    return NUM_PENTAGONS
end

"""
    getPentagons(res::Int) -> Tuple{H3Error, Vector{H3Index}}

Get all pentagon cell indexes at the given resolution.

See also the H3 C API: [`getPentagons`](https://h3geo.org/docs/api/misc#getpentagons)
"""
function getPentagons(res::Int)::Tuple{H3Error, Vector{H3Index}}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, H3Index[])
    end
    out = H3Index[]
    for bc in 0:(NUM_BASE_CELLS - 1)
        if _isBaseCellPentagon(bc)
            push!(out, setH3Index(res, bc, CENTER_DIGIT))
        end
    end
    return (E_SUCCESS, out)
end

function res0CellCount_h3()::Int
    return NUM_BASE_CELLS
end

"""
    getRes0Cells() -> Tuple{H3Error, Vector{H3Index}}

Get all resolution 0 (base cell) indexes.

See also the H3 C API: [`getRes0Cells`](https://h3geo.org/docs/api/misc#getres0cells)
"""
function getRes0Cells()::Tuple{H3Error, Vector{H3Index}}
    out = Vector{H3Index}(undef, NUM_BASE_CELLS)
    for bc in 0:(NUM_BASE_CELLS - 1)
        h = H3_INIT
        h = h3_set_mode(h, H3_CELL_MODE)
        h = h3_set_base_cell(h, bc)
        out[bc + 1] = h
    end
    return (E_SUCCESS, out)
end

"""
    stringToH3(str::AbstractString) -> Tuple{H3Error, H3Index}

Parse a hexadecimal string to an H3 index.

See also the H3 C API: [`stringToH3`](https://h3geo.org/docs/api/inspection#stringtoh3)
"""
function stringToH3(str::AbstractString)::Tuple{H3Error, H3Index}
    try
        h = parse(UInt64, str, base=16)
        return (E_SUCCESS, h)
    catch
        return (E_FAILED, H3_NULL)
    end
end

"""
    h3ToString(h::H3Index) -> String

Convert an H3 index to its hexadecimal string representation.

See also the H3 C API: [`h3ToString`](https://h3geo.org/docs/api/inspection#h3tostring)
"""
function h3ToString(h::H3Index)::String
    return string(h, base=16)
end

"""
    constructCell(res::Int, baseCellNumber::Int, digits::Vector{Int}) -> Tuple{H3Error, H3Index}

Construct an H3 cell index from resolution, base cell number, and digit array.
"""
function constructCell(res::Int, baseCellNumber::Int, digits::Vector{Int})::Tuple{H3Error, H3Index}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, H3_NULL)
    end
    if baseCellNumber < 0 || baseCellNumber >= NUM_BASE_CELLS
        return (E_BASE_CELL_DOMAIN, H3_NULL)
    end
    h = H3_INIT
    h = h3_set_mode(h, H3_CELL_MODE)
    h = h3_set_resolution(h, res)
    h = h3_set_base_cell(h, baseCellNumber)

    isPent = _isBaseCellPentagon(baseCellNumber)
    for r in 1:res
        d = digits[r]
        if d < 0 || d >= 7
            return (E_DIGIT_DOMAIN, H3_NULL)
        end
        if isPent
            if d == 0
                # still on pentagon center
            elseif d == 1
                return (E_DELETED_DIGIT, H3_NULL)
            else
                isPent = false
            end
        end
        h = h3_set_index_digit(h, r, d)
    end
    return (E_SUCCESS, h)
end

"""
    getIcosahedronFaces(h::H3Index) -> Tuple{H3Error, Vector{Int}}

Get all icosahedron faces intersected by a cell.

See also the H3 C API: [`getIcosahedronFaces`](https://h3geo.org/docs/api/inspection#geticosahedronfaces)
"""
function getIcosahedronFaces(h::H3Index)::Tuple{H3Error, Vector{Int}}
    res = h3_get_resolution(h)
    isPent = isPentagon(h)

    if isPent && !isResolutionClassIII(res)
        childPentagon = makeDirectChild(h, 0)
        return getIcosahedronFaces(childPentagon)
    end

    err, fijk = _h3ToFaceIjk(h)
    if err != E_SUCCESS
        return (err, Int[])
    end

    err2, faceCount = maxFaceCount(h)
    if err2 != E_SUCCESS
        return (err2, Int[])
    end

    if isPent
        vertexCount = NUM_PENT_VERTS
        fijkVerts, adjRes = _faceIjkPentToVerts(fijk, res)
    else
        vertexCount = NUM_HEX_VERTS
        fijkVerts, adjRes = _faceIjkToVerts(fijk, res)
    end

    out = fill(INVALID_FACE, faceCount)

    for i in 1:vertexCount
        vert = fijkVerts[i]
        if isPent
            _, vert = _adjustPentVertOverage(vert, adjRes)
        else
            _, vert = _adjustOverageClassII(vert, adjRes, 0, 1)
        end

        face = Int(vert.face)
        pos = 1
        while pos <= faceCount && out[pos] != INVALID_FACE && out[pos] != face
            pos += 1
        end
        if pos > faceCount
            return (E_FAILED, Int[])
        end
        out[pos] = face
    end
    return (E_SUCCESS, out)
end

# Compact/uncompact operations
"""
    compactCells(h3Set::Vector{H3Index}) -> Tuple{H3Error, Vector{H3Index}}

Compact a set of cells into the smallest equivalent set.

See also the H3 C API: [`compactCells`](https://h3geo.org/docs/api/hierarchy#compactcells)
"""
function compactCells(h3Set::Vector{H3Index})::Tuple{H3Error, Vector{H3Index}}
    numHexes = length(h3Set)
    if numHexes == 0
        return (E_SUCCESS, H3Index[])
    end

    res = h3_get_resolution(h3Set[1])
    if res == 0
        return (E_SUCCESS, copy(h3Set))
    end

    remainingHexes = copy(h3Set)
    compactedSet = H3Index[]

    while !isempty(remainingHexes)
        res = h3_get_resolution(remainingHexes[1])
        parentRes = res - 1

        if parentRes < 0
            append!(compactedSet, remainingHexes)
            break
        end

        parentCounts = Dict{H3Index, Int}()
        for idx in remainingHexes
            if idx != H3_NULL
                if h3_get_reserved_bits(idx) != 0
                    return (E_CELL_INVALID, H3Index[])
                end
                err, parent = cellToParent(idx, parentRes)
                if err != E_SUCCESS
                    return (err, H3Index[])
                end
                parentCounts[parent] = get(parentCounts, parent, 0) + 1
            end
        end

        compactableParents = Set{H3Index}()
        for (parent, count) in parentCounts
            limitCount = isPentagon(parent) ? 6 : 7
            if count > limitCount
                return (E_DUPLICATE_INPUT, H3Index[])
            end
            needed = isPentagon(parent) ? 6 : 7
            if count == needed
                push!(compactableParents, parent)
            end
        end

        uncompactable = H3Index[]
        nextRemaining = H3Index[]

        for idx in remainingHexes
            if idx != H3_NULL
                if parentRes >= 0
                    err, parent = cellToParent(idx, parentRes)
                    if err != E_SUCCESS
                        return (err, H3Index[])
                    end
                    if parent in compactableParents
                        continue
                    end
                end
                push!(uncompactable, idx)
            end
        end

        append!(compactedSet, uncompactable)
        remainingHexes = collect(compactableParents)

        if isempty(remainingHexes)
            break
        end
    end

    return (E_SUCCESS, compactedSet)
end

"""
    uncompactCellsSize(compactedSet::Vector{H3Index}, res::Int) -> Tuple{H3Error, Int64}

Get the number of cells that would result from uncompacting to the given resolution.

See also the H3 C API: [`uncompactCellsSize`](https://h3geo.org/docs/api/hierarchy#uncompactcellssize)
"""
function uncompactCellsSize(compactedSet::Vector{H3Index}, res::Int)::Tuple{H3Error, Int64}
    numOut = Int64(0)
    for h in compactedSet
        if h == H3_NULL
            continue
        end
        err, childSize = cellToChildrenSize(h, res)
        if err != E_SUCCESS
            return (E_RES_MISMATCH, Int64(0))
        end
        numOut += childSize
    end
    return (E_SUCCESS, numOut)
end

"""
    uncompactCells(compactedSet::Vector{H3Index}, res::Int) -> Tuple{H3Error, Vector{H3Index}}

Uncompact a set of cells to the given resolution.

See also the H3 C API: [`uncompactCells`](https://h3geo.org/docs/api/hierarchy#uncompactcells)
"""
function uncompactCells(compactedSet::Vector{H3Index}, res::Int)::Tuple{H3Error, Vector{H3Index}}
    err, totalSize = uncompactCellsSize(compactedSet, res)
    if err != E_SUCCESS
        return (err, H3Index[])
    end
    out = Vector{H3Index}(undef, totalSize)
    idx = 1
    for h in compactedSet
        if h == H3_NULL
            continue
        end
        hRes = h3_get_resolution(h)
        if res < hRes
            return (E_RES_MISMATCH, H3Index[])
        end
        _uncompactCellRecursive!(h, res, out, idx)
        err2, childSize = cellToChildrenSize(h, res)
        idx += childSize
    end
    return (E_SUCCESS, out)
end

function _uncompactCellRecursive!(h::H3Index, targetRes::Int, out::Vector{H3Index}, startIdx::Int)
    hRes = h3_get_resolution(h)
    if hRes == targetRes
        out[startIdx] = h
        return
    end
    isPent = isPentagon(h)
    idx = startIdx
    for d in 0:6
        if isPent && d == 1
            continue
        end
        child = makeDirectChild(h, d)
        if h3_get_resolution(child) == targetRes
            out[idx] = child
            idx += 1
        else
            err, childSize = cellToChildrenSize(child, targetRes)
            _uncompactCellRecursive!(child, targetRes, out, idx)
            idx += childSize
        end
    end
end

"""
    cellToChildPos(child::H3Index, parentRes::Int) -> Tuple{H3Error, Int64}

Get the position of a child cell within its parent's ordered children.

See also the H3 C API: [`cellToChildPos`](https://h3geo.org/docs/api/hierarchy#celltochildpos)
"""
function cellToChildPos(child::H3Index, parentRes::Int)::Tuple{H3Error, Int64}
    childRes = h3_get_resolution(child)
    err, originalParent = cellToParent(child, parentRes)
    if err != E_SUCCESS
        return (err, Int64(0))
    end

    parent = originalParent
    parentIsPentagon = isPentagon(parent)

    pos = Int64(0)
    if parentIsPentagon
        for res in childRes:-1:(parentRes + 1)
            err2, parent = cellToParent(child, res - 1)
            if err2 != E_SUCCESS
                return (err2, Int64(0))
            end
            parentIsPentagon = isPentagon(parent)
            rawDigit = Int(h3_get_index_digit(child, res))
            if rawDigit == Int(INVALID_DIGIT) ||
               (parentIsPentagon && rawDigit == Int(K_AXES_DIGIT))
                return (E_CELL_INVALID, Int64(0))
            end
            digit = parentIsPentagon && rawDigit > 0 ? rawDigit - 1 : rawDigit
            if digit != 0
                hexChildCount = _ipow(Int64(7), Int64(childRes - res))
                pentOffset = parentIsPentagon ?
                    Int64(1) + (Int64(5) * (hexChildCount - Int64(1))) ÷ Int64(6) :
                    hexChildCount
                pos += pentOffset + Int64(digit - 1) * hexChildCount
            end
        end
    else
        for res in childRes:-1:(parentRes + 1)
            digit = Int(h3_get_index_digit(child, res))
            if digit == Int(INVALID_DIGIT)
                return (E_CELL_INVALID, Int64(0))
            end
            pos += Int64(digit) * _ipow(Int64(7), Int64(childRes - res))
        end
    end

    return (E_SUCCESS, pos)
end

"""
    childPosToCell(childPos::Int64, parent::H3Index, childRes::Int) -> Tuple{H3Error, H3Index}

Convert a child position to the corresponding child cell.

See also the H3 C API: [`childPosToCell`](https://h3geo.org/docs/api/hierarchy#childpostocell)
"""
function childPosToCell(childPos::Int64, parent::H3Index, childRes::Int)::Tuple{H3Error, H3Index}
    if childRes < 0 || childRes > MAX_H3_RES
        return (E_RES_DOMAIN, H3_NULL)
    end
    parentRes = h3_get_resolution(parent)
    if childRes < parentRes
        return (E_RES_MISMATCH, H3_NULL)
    end

    err, maxChildCount = cellToChildrenSize(parent, childRes)
    if err != E_SUCCESS
        return (err, H3_NULL)
    end
    if childPos < 0 || childPos >= maxChildCount
        return (E_DOMAIN, H3_NULL)
    end

    resOffset = childRes - parentRes
    child = h3_set_resolution(parent, childRes)
    idx = childPos

    if isPentagon(parent)
        inPent = true
        for res in 1:resOffset
            resWidth = _ipow(Int64(7), Int64(resOffset - res))
            if inPent
                pentWidth = Int64(1) + (Int64(5) * (resWidth - Int64(1))) ÷ Int64(6)
                if idx < pentWidth
                    child = h3_set_index_digit(child, parentRes + res, 0)
                else
                    idx -= pentWidth
                    inPent = false
                    child = h3_set_index_digit(child, parentRes + res, Int(idx ÷ resWidth) + 2)
                    idx = idx % resWidth
                end
            else
                child = h3_set_index_digit(child, parentRes + res, Int(idx ÷ resWidth))
                idx = idx % resWidth
            end
        end
    else
        for res in 1:resOffset
            resWidth = _ipow(Int64(7), Int64(resOffset - res))
            child = h3_set_index_digit(child, parentRes + res, Int(idx ÷ resWidth))
            idx = idx % resWidth
        end
    end

    return (E_SUCCESS, child)
end

"""
    cellToChildren(h::H3Index, childRes::Int) -> Tuple{H3Error, Vector{H3Index}}

Get all children of a cell at a given finer resolution.

See also the H3 C API: [`cellToChildren`](https://h3geo.org/docs/api/hierarchy#celltochildren)
"""
function cellToChildren(h::H3Index, childRes::Int)::Tuple{H3Error, Vector{H3Index}}
    err, numChildren = cellToChildrenSize(h, childRes)
    if err != E_SUCCESS
        return (err, H3Index[])
    end
    out = Vector{H3Index}(undef, numChildren)
    for i in Int64(0):(numChildren - Int64(1))
        err2, child = childPosToCell(i, h, childRes)
        if err2 != E_SUCCESS
            return (err2, H3Index[])
        end
        out[i + 1] = child
    end
    return (E_SUCCESS, out)
end

"""
    isValidIndex(h::H3Index) -> Bool

Check whether an H3 index is valid (cell, directed edge, or vertex).

See also the H3 C API: [`isValidIndex`](https://h3geo.org/docs/api/inspection#isvalidindex)
"""
function isValidIndex(h::H3Index)::Bool
    return isValidCell(h) || isValidDirectedEdge(h) || isValidVertex(h)
end
