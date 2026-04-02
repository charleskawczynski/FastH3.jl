# Polygon utility functions

function maxPolygonToCellsSize(polygon::GeoPolygon, res::Int, flags::UInt32)::Tuple{H3Error, Int64}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, Int64(0))
    end
    bbox = bboxFromGeoLoop(polygon.geoloop)
    err, numHexes = _bboxHexEstimate(bbox, res)
    if err != E_SUCCESS
        return (err, Int64(0))
    end
    return (E_SUCCESS, numHexes)
end

function _bboxHexEstimate(bbox::BBox, res::Int)::Tuple{H3Error, Int64}
    err, pentagons = getPentagons(res)
    if err != E_SUCCESS
        return (err, Int64(0))
    end
    p1 = LatLng(bbox.north, bbox.east)
    p2 = LatLng(bbox.south, bbox.west)
    d = greatCircleDistanceKm(p1, p2)
    lngDiff = abs(p1.lng - p2.lng)
    latDiff = abs(p1.lat - p2.lat)
    if lngDiff == 0.0 || latDiff == 0.0
        return (E_FAILED, Int64(0))
    end
    
    err2, pentArea = getHexagonAreaAvgKm2(res)
    if err2 != E_SUCCESS
        return (err2, Int64(0))
    end
    
    length_val = max(lngDiff, latDiff)
    width_val = min(lngDiff, latDiff)
    ratio = length_val / width_val
    a = d * d / min(3.0, ratio)
    estimateDouble = ceil(a / (pentArea * 0.8))
    if !isfinite(estimateDouble)
        return (E_FAILED, Int64(0))
    end
    estimate = max(Int64(1), Int64(trunc(estimateDouble)))
    return (E_SUCCESS, estimate)
end

function polygonToCells(polygon::GeoPolygon, res::Int, flags::UInt32)::Tuple{H3Error, Vector{H3Index}}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, H3Index[])
    end
    
    bbox = bboxFromGeoLoop(polygon.geoloop)
    
    err, maxSize = maxPolygonToCellsSize(polygon, res, flags)
    if err != E_SUCCESS
        return (err, H3Index[])
    end
    
    out = H3Index[]
    
    err2, minCell = latLngToCell(LatLng(bbox.south, bbox.west), res)
    if err2 != E_SUCCESS
        return (err2, H3Index[])
    end
    
    err3, maxCell = latLngToCell(LatLng(bbox.north, bbox.east), res)
    if err3 != E_SUCCESS
        return (err3, H3Index[])
    end
    
    err4, numCells_estimate = maxPolygonToCellsSize(polygon, res, flags)
    if err4 != E_SUCCESS
        return (err4, H3Index[])
    end
    
    # BFS from center cell, only adding cells whose centers are in the polygon
    center = bboxCenter(bbox)
    err5, centerCell = latLngToCell(center, res)
    if err5 != E_SUCCESS
        return (err5, H3Index[])
    end
    
    visited = Set{H3Index}()
    queue = H3Index[centerCell]
    push!(visited, centerCell)
    
    while !isempty(queue)
        current = popfirst!(queue)
        err6, center_ll = cellToLatLng(current)
        if err6 != E_SUCCESS
            continue
        end
        
        if pointInsideGeoLoop(polygon.geoloop, bbox, center_ll)
            inHole = false
            for hi in 1:polygon.numHoles
                holeBbox = bboxFromGeoLoop(polygon.holes[hi])
                if pointInsideGeoLoop(polygon.holes[hi], holeBbox, center_ll)
                    inHole = true
                    break
                end
            end
            if !inHole
                push!(out, current)
                err7, neighbors = gridDisk(current, 1)
                if err7 == E_SUCCESS
                    for n in neighbors
                        if n != H3_NULL && !(n in visited)
                            push!(visited, n)
                            push!(queue, n)
                        end
                    end
                end
            end
        end
        
        if length(visited) > numCells_estimate * 2
            break
        end
    end
    
    return (E_SUCCESS, out)
end
