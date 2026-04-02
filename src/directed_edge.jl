# Directed edge functions
# Translated from H3 library directedEdge.c

"""
Check whether an H3 index represents a valid directed edge.
"""
function isValidDirectedEdge(edge::H3Index)::Bool
    if h3_get_mode(edge) != H3_DIRECTEDEDGE_MODE
        return false
    end
    neighborDirection = h3_get_reserved_bits(edge)
    if neighborDirection < Int(K_AXES_DIGIT) || neighborDirection > Int(IJ_AXES_DIGIT)
        return false
    end
    origin = h3_set_mode(edge, H3_CELL_MODE)
    origin = h3_set_reserved_bits(origin, 0)
    return isValidCell(origin)
end

"""
Get the origin cell from a directed edge.
"""
function getDirectedEdgeOrigin(edge::H3Index)::Tuple{H3Error, H3Index}
    if h3_get_mode(edge) != H3_DIRECTEDEDGE_MODE
        return (E_DIR_EDGE_INVALID, H3_NULL)
    end
    origin = h3_set_mode(edge, H3_CELL_MODE)
    origin = h3_set_reserved_bits(origin, 0)
    return (E_SUCCESS, origin)
end

"""
Get the destination cell from a directed edge by finding the neighbor
of the origin in the encoded direction.
"""
function getDirectedEdgeDestination(edge::H3Index)::Tuple{H3Error, H3Index}
    if h3_get_mode(edge) != H3_DIRECTEDEDGE_MODE
        return (E_DIR_EDGE_INVALID, H3_NULL)
    end
    direction = h3_get_reserved_bits(edge)
    if direction < Int(K_AXES_DIGIT) || direction > Int(IJ_AXES_DIGIT)
        return (E_DIR_EDGE_INVALID, H3_NULL)
    end

    origin = h3_set_mode(edge, H3_CELL_MODE)
    origin = h3_set_reserved_bits(origin, 0)

    err, dest, _ = _h3NeighborRotations(origin, Direction(direction), 0)
    if err != E_SUCCESS
        return (err, H3_NULL)
    end
    return (E_SUCCESS, dest)
end

"""
Get both the origin and destination cells from a directed edge.
"""
function directedEdgeToCells(edge::H3Index)::Tuple{H3Error, H3Index, H3Index}
    err1, origin = getDirectedEdgeOrigin(edge)
    if err1 != E_SUCCESS
        return (err1, H3_NULL, H3_NULL)
    end
    err2, dest = getDirectedEdgeDestination(edge)
    if err2 != E_SUCCESS
        return (err2, H3_NULL, H3_NULL)
    end
    return (E_SUCCESS, origin, dest)
end

"""
Determine whether two H3 cells are neighbors (share an edge).
"""
function areNeighborCells(origin::H3Index, destination::H3Index)::Tuple{H3Error, Bool}
    if h3_get_mode(origin) != H3_CELL_MODE
        return (E_CELL_INVALID, false)
    end
    if h3_get_mode(destination) != H3_CELL_MODE
        return (E_CELL_INVALID, false)
    end
    if origin == destination
        return (E_SUCCESS, false)
    end

    originRes = h3_get_resolution(origin)
    destRes = h3_get_resolution(destination)
    if originRes != destRes
        return (E_RES_MISMATCH, false)
    end

    err, ij = cellToLocalIj(origin, destination, UInt32(0))
    if err != E_SUCCESS
        err2, ij2 = cellToLocalIj(destination, origin, UInt32(0))
        if err2 != E_SUCCESS
            return (err2, false)
        end
        ij = ij2
    end

    err3, ijk = ijToIjk(ij)
    if err3 != E_SUCCESS
        return (err3, false)
    end

    return (E_SUCCESS, ijkDistance(ijk, CoordIJK(0, 0, 0)) == 1)
end

"""
Create a directed edge H3 index from an origin cell to a destination cell.
The cells must be neighbors.
"""
function cellsToDirectedEdge(origin::H3Index, destination::H3Index)::Tuple{H3Error, H3Index}
    err, isNeighbor = areNeighborCells(origin, destination)
    if err != E_SUCCESS
        return (err, H3_NULL)
    end
    if !isNeighbor
        return (E_NOT_NEIGHBORS, H3_NULL)
    end

    edge = h3_set_mode(origin, H3_DIRECTEDEDGE_MODE)

    isPent = isPentagon(origin)
    for dir_i in Int(K_AXES_DIGIT):Int(IJ_AXES_DIGIT)
        if isPent && dir_i == Int(K_AXES_DIGIT)
            continue
        end
        d = Direction(dir_i)
        nerr, neighbor, _ = _h3NeighborRotations(origin, d, 0)
        if nerr == E_SUCCESS && neighbor == destination
            edge = h3_set_reserved_bits(edge, dir_i)
            return (E_SUCCESS, edge)
        end
    end

    return (E_NOT_NEIGHBORS, H3_NULL)
end

"""
Get all directed edges originating from a cell.
Returns 6 edges for hexagons, 5 for pentagons.
"""
function originToDirectedEdges(origin::H3Index)::Tuple{H3Error, Vector{H3Index}}
    if !isValidCell(origin)
        return (E_CELL_INVALID, H3Index[])
    end

    isPent = isPentagon(origin)
    edges = H3Index[]

    for dir_i in Int(K_AXES_DIGIT):Int(IJ_AXES_DIGIT)
        if isPent && dir_i == Int(K_AXES_DIGIT)
            continue
        end
        edge = h3_set_mode(origin, H3_DIRECTEDEDGE_MODE)
        edge = h3_set_reserved_bits(edge, dir_i)
        push!(edges, edge)
    end

    return (E_SUCCESS, edges)
end

"""
Get the boundary of a directed edge (the vertices shared between origin
and destination cells).
"""
function directedEdgeToBoundary(edge::H3Index)::Tuple{H3Error, CellBoundary}
    err, origin, dest = directedEdgeToCells(edge)
    if err != E_SUCCESS
        return (err, CellBoundary())
    end

    err1, cb1 = cellToBoundary(origin)
    if err1 != E_SUCCESS
        return (err1, CellBoundary())
    end

    err2, cb2 = cellToBoundary(dest)
    if err2 != E_SUCCESS
        return (err2, CellBoundary())
    end

    shared = LatLng[]
    for i in 1:cb1.numVerts
        for j in 1:cb2.numVerts
            v1 = cb1.verts[i]
            v2 = cb2.verts[j]
            if geoAlmostEqualThreshold(v1, v2, EPSILON_RAD * 100.0)
                push!(shared, v1)
                break
            end
        end
    end

    numVerts = min(Int32(length(shared)), Int32(MAX_CELL_BNDRY_VERTS))
    verts_tup = ntuple(i -> i <= length(shared) ? shared[i] : LatLng(), MAX_CELL_BNDRY_VERTS)

    return (E_SUCCESS, CellBoundary(numVerts, verts_tup))
end

"""
Get the reversed directed edge (swapping origin and destination).
"""
function reverseDirectedEdge(edge::H3Index)::Tuple{H3Error, H3Index}
    if !isValidDirectedEdge(edge)
        return (E_DIR_EDGE_INVALID, H3_NULL)
    end

    err1, origin = getDirectedEdgeOrigin(edge)
    if err1 != E_SUCCESS
        return (err1, H3_NULL)
    end

    err2, dest = getDirectedEdgeDestination(edge)
    if err2 != E_SUCCESS
        return (err2, H3_NULL)
    end

    return cellsToDirectedEdge(dest, origin)
end
