# Vertex functions
# Translated from H3 library vertex.c

struct PentagonDirectionFaces
    baseCell::Int32
    faces::NTuple{5, Int32}
end

const pentagonDirectionFaces = [
    PentagonDirectionFaces(4,   (4, 0, 2, 1, 3)),
    PentagonDirectionFaces(14,  (6, 11, 2, 7, 1)),
    PentagonDirectionFaces(24,  (5, 10, 1, 6, 0)),
    PentagonDirectionFaces(38,  (7, 12, 3, 8, 2)),
    PentagonDirectionFaces(49,  (9, 14, 0, 5, 4)),
    PentagonDirectionFaces(58,  (8, 13, 4, 9, 3)),
    PentagonDirectionFaces(63,  (11, 6, 15, 10, 16)),
    PentagonDirectionFaces(72,  (12, 7, 16, 11, 17)),
    PentagonDirectionFaces(83,  (10, 5, 19, 14, 15)),
    PentagonDirectionFaces(97,  (13, 8, 17, 12, 18)),
    PentagonDirectionFaces(107, (14, 9, 18, 13, 19)),
    PentagonDirectionFaces(117, (15, 19, 17, 18, 16)),
]

# Direction-to-vertex mappings for hexagons and pentagons.
# For hexagons, vertex N is between edges in directions
# directionForVertexNum[N+1] and the next direction CCW.
const _directionForVertexNumHex = Direction[
    J_AXES_DIGIT, JK_AXES_DIGIT, K_AXES_DIGIT,
    IK_AXES_DIGIT, I_AXES_DIGIT, IJ_AXES_DIGIT,
]

const _directionForVertexNumPent = Direction[
    J_AXES_DIGIT, JK_AXES_DIGIT,
    IK_AXES_DIGIT, I_AXES_DIGIT, IJ_AXES_DIGIT,
]

"""
    isValidVertex(vertex::H3Index) -> Bool

Check whether an H3 index represents a valid vertex.

See also the H3 C API: [`isValidVertex`](https://h3geo.org/docs/api/vertex#isvalidvertex)
"""
function isValidVertex(vertex::H3Index)::Bool
    if h3_get_mode(vertex) != H3_VERTEX_MODE
        return false
    end
    vertexNum = h3_get_reserved_bits(vertex)
    owner = h3_set_mode(vertex, H3_CELL_MODE)
    owner = h3_set_reserved_bits(owner, 0)
    if !isValidCell(owner)
        return false
    end
    numVerts = isPentagon(owner) ? NUM_PENT_VERTS : NUM_HEX_VERTS
    if vertexNum < 0 || vertexNum >= numVerts
        return false
    end
    return true
end

"""
    cellToVertex(origin::H3Index, vertexNum::Int) -> (H3Error, H3Index)

Get a vertex H3 index for vertex number `vertexNum` (0-indexed) of cell `origin`.

See also the H3 C API: [`cellToVertex`](https://h3geo.org/docs/api/vertex#celltovertex)
"""
function cellToVertex(origin::H3Index, vertexNum::Int)::Tuple{H3Error, H3Index}
    if !isValidCell(origin)
        return (E_CELL_INVALID, H3_NULL)
    end

    numVerts = isPentagon(origin) ? NUM_PENT_VERTS : NUM_HEX_VERTS
    if vertexNum < 0 || vertexNum >= numVerts
        return (E_DOMAIN, H3_NULL)
    end

    vertex = h3_set_mode(origin, H3_VERTEX_MODE)
    vertex = h3_set_reserved_bits(vertex, vertexNum)

    return (E_SUCCESS, vertex)
end

"""
    cellToVertexes(origin::H3Index) -> (H3Error, Vector{H3Index})

Get all vertex H3 indexes for a cell (5 for pentagons, 6 for hexagons).

See also the H3 C API: [`cellToVertexes`](https://h3geo.org/docs/api/vertex#celltovertexes)
"""
function cellToVertexes(origin::H3Index)::Tuple{H3Error, Vector{H3Index}}
    if !isValidCell(origin)
        return (E_CELL_INVALID, H3Index[])
    end

    numVerts = isPentagon(origin) ? NUM_PENT_VERTS : NUM_HEX_VERTS
    verts = Vector{H3Index}(undef, numVerts)

    for v in 0:(numVerts - 1)
        err, vert = cellToVertex(origin, v)
        if err != E_SUCCESS
            return (err, H3Index[])
        end
        verts[v + 1] = vert
    end

    return (E_SUCCESS, verts)
end

"""
    vertexToLatLng(vertex::H3Index) -> (H3Error, LatLng)

Get the latitude/longitude of a vertex H3 index.

See also the H3 C API: [`vertexToLatLng`](https://h3geo.org/docs/api/vertex#vertextolatlng)
"""
function vertexToLatLng(vertex::H3Index)::Tuple{H3Error, LatLng}
    if !isValidVertex(vertex)
        return (E_VERTEX_INVALID, LatLng())
    end

    vertexNum = h3_get_reserved_bits(vertex)
    owner = h3_set_mode(vertex, H3_CELL_MODE)
    owner = h3_set_reserved_bits(owner, 0)

    err, cb = cellToBoundary(owner)
    if err != E_SUCCESS
        return (err, LatLng())
    end

    if vertexNum >= cb.numVerts
        return (E_VERTEX_INVALID, LatLng())
    end

    return (E_SUCCESS, cb.verts[vertexNum + 1])
end
