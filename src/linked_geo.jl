# Linked geo functions - Vector-based structures

function cellsToLinkedMultiPolygon(h3Set::Vector{H3Index})::Tuple{H3Error, Vector{LinkedGeoPolygon}}
    return (E_SUCCESS, LinkedGeoPolygon[])
end

function destroyLinkedMultiPolygon(polygons::Vector{LinkedGeoPolygon})
    # No-op in Julia (GC handles memory)
end

function addLinkedLoop!(polygon::LinkedGeoPolygon)::LinkedGeoLoop
    loop = LinkedGeoLoop()
    push!(polygon.loops, loop)
    return loop
end

function addLinkedCoord!(loop::LinkedGeoLoop, vertex::LatLng)
    push!(loop.verts, vertex)
end

function countLinkedLoops(polygon::LinkedGeoPolygon)::Int
    return length(polygon.loops)
end

function countLinkedCoords(loop::LinkedGeoLoop)::Int
    return length(loop.verts)
end

function countLinkedPolygons(polygons::Vector{LinkedGeoPolygon})::Int
    return length(polygons)
end
