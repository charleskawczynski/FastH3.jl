# Geographic bounding box functions

struct BBox
    north::Float64
    south::Float64
    east::Float64
    west::Float64
end
BBox() = BBox(0.0, 0.0, 0.0, 0.0)

function bboxWidthRads(bbox::BBox)::Float64
    return bboxIsTransmeridian(bbox) ? bbox.east - bbox.west + M_2PI :
           bbox.east - bbox.west
end

function bboxHeightRads(bbox::BBox)::Float64
    return bbox.north - bbox.south
end

function bboxIsTransmeridian(bbox::BBox)::Bool
    return bbox.east < bbox.west
end

function bboxCenter(bbox::BBox)::LatLng
    lat = (bbox.north + bbox.south) * 0.5
    east = bboxIsTransmeridian(bbox) ? bbox.east + M_2PI : bbox.east
    lng = constrainLng((east + bbox.west) * 0.5)
    return LatLng(lat, lng)
end

function bboxContains(bbox::BBox, point::LatLng)::Bool
    return point.lat >= bbox.south && point.lat <= bbox.north &&
           (bboxIsTransmeridian(bbox) ?
                (point.lng >= bbox.west || point.lng <= bbox.east) :
                (point.lng >= bbox.west && point.lng <= bbox.east))
end

function bboxOverlapsBBox(a::BBox, b::BBox)::Bool
    if a.north < b.south || a.south > b.north
        return false
    end
    aNorm, bNorm = bboxNormalization(a, b)
    if normalizeLng(a.east, aNorm) < normalizeLng(b.west, bNorm) ||
       normalizeLng(a.west, aNorm) > normalizeLng(b.east, bNorm)
        return false
    end
    return true
end

function bboxContainsBBox(a::BBox, b::BBox)::Bool
    if a.north < b.north || a.south > b.south
        return false
    end
    aNorm, bNorm = bboxNormalization(a, b)
    return normalizeLng(a.west, aNorm) <= normalizeLng(b.west, bNorm) &&
           normalizeLng(a.east, aNorm) >= normalizeLng(b.east, bNorm)
end

function bboxEquals(b1::BBox, b2::BBox)::Bool
    return b1.north == b2.north && b1.south == b2.south &&
           b1.east == b2.east && b1.west == b2.west
end

function bboxToCellBoundary(bbox::BBox)::CellBoundary
    verts = ntuple(i -> begin
        if i == 1
            LatLng(bbox.north, bbox.east)
        elseif i == 2
            LatLng(bbox.north, bbox.west)
        elseif i == 3
            LatLng(bbox.south, bbox.west)
        elseif i == 4
            LatLng(bbox.south, bbox.east)
        else
            LatLng()
        end
    end, MAX_CELL_BNDRY_VERTS)
    return CellBoundary(Int32(4), verts)
end

function scaleBBox(bbox::BBox, scale::Float64)::BBox
    width = bboxWidthRads(bbox)
    height = bboxHeightRads(bbox)
    widthBuffer = (width * scale - width) * 0.5
    heightBuffer = (height * scale - height) * 0.5
    north = bbox.north + heightBuffer
    if north > M_PI_2
        north = M_PI_2
    end
    south = bbox.south - heightBuffer
    if south < -M_PI_2
        south = -M_PI_2
    end
    east = bbox.east + widthBuffer
    if east > M_PI
        east -= M_2PI
    end
    if east < -M_PI
        east += M_2PI
    end
    west = bbox.west - widthBuffer
    if west > M_PI
        west -= M_2PI
    end
    if west < -M_PI
        west += M_2PI
    end
    return BBox(north, south, east, west)
end

function bboxNormalization(a::BBox, b::BBox)::Tuple{LongitudeNormalization, LongitudeNormalization}
    aIsTransmeridian = bboxIsTransmeridian(a)
    bIsTransmeridian = bboxIsTransmeridian(b)
    aToBTrendsEast = a.west - b.east < b.west - a.east

    aNormalization = if !aIsTransmeridian
        NORMALIZE_NONE
    elseif bIsTransmeridian
        NORMALIZE_EAST
    elseif aToBTrendsEast
        NORMALIZE_EAST
    else
        NORMALIZE_WEST
    end

    bNormalization = if !bIsTransmeridian
        NORMALIZE_NONE
    elseif aIsTransmeridian
        NORMALIZE_EAST
    elseif aToBTrendsEast
        NORMALIZE_WEST
    else
        NORMALIZE_EAST
    end

    return (aNormalization, bNormalization)
end

function bboxFromGeoLoop(loop::GeoLoop)::BBox
    if loop.numVerts == 0
        return BBox(0.0, 0.0, 0.0, 0.0)
    end

    north = -Inf
    south = Inf
    east = -Inf
    west = Inf
    minPosLng = Inf
    maxNegLng = -Inf
    isTransmeridian = false

    for idx in 1:loop.numVerts
        coord = loop.verts[idx]
        nextIdx = mod1(idx + 1, loop.numVerts)
        next = loop.verts[nextIdx]

        if coord.lat < south
            south = coord.lat
        end
        if coord.lng < west
            west = coord.lng
        end
        if coord.lat > north
            north = coord.lat
        end
        if coord.lng > east
            east = coord.lng
        end
        if coord.lng > 0 && coord.lng < minPosLng
            minPosLng = coord.lng
        end
        if coord.lng < 0 && coord.lng > maxNegLng
            maxNegLng = coord.lng
        end
        if abs(coord.lng - next.lng) > M_PI
            isTransmeridian = true
        end
    end

    if isTransmeridian
        east = maxNegLng
        west = minPosLng
    end

    return BBox(north, south, east, west)
end

function pointInsideGeoLoop(loop::GeoLoop, bbox::BBox, coord::LatLng)::Bool
    if !bboxContains(bbox, coord)
        return false
    end
    isTransmeridian = bboxIsTransmeridian(bbox)
    contains = false

    lat = coord.lat
    lng = isTransmeridian && coord.lng < 0 ? coord.lng + M_2PI : coord.lng

    for idx in 1:loop.numVerts
        a = loop.verts[idx]
        nextIdx = mod1(idx + 1, loop.numVerts)
        b = loop.verts[nextIdx]

        if a.lat > b.lat
            a, b = b, a
        end

        if lat == a.lat || lat == b.lat
            lat += eps(Float64)
        end

        if lat < a.lat || lat > b.lat
            continue
        end

        aLng = isTransmeridian && a.lng < 0 ? a.lng + M_2PI : a.lng
        bLng = isTransmeridian && b.lng < 0 ? b.lng + M_2PI : b.lng

        if aLng == lng || bLng == lng
            lng -= eps(Float64)
        end

        ratio = (lat - a.lat) / (b.lat - a.lat)
        testLng_val = aLng + (bLng - aLng) * ratio
        if isTransmeridian && testLng_val < 0
            testLng_val += M_2PI
        end

        if testLng_val > lng
            contains = !contains
        end
    end

    return contains
end

function isClockwiseGeoLoop(loop::GeoLoop)::Bool
    return _isClockwiseNormalizedGeoLoop(loop, false)
end

function _isClockwiseNormalizedGeoLoop(loop::GeoLoop, isTransmeridian::Bool)::Bool
    sum = 0.0
    for idx in 1:loop.numVerts
        a = loop.verts[idx]
        nextIdx = mod1(idx + 1, loop.numVerts)
        b = loop.verts[nextIdx]

        if !isTransmeridian && abs(a.lng - b.lng) > M_PI
            return _isClockwiseNormalizedGeoLoop(loop, true)
        end

        aLng = isTransmeridian && a.lng < 0 ? a.lng + M_2PI : a.lng
        bLng = isTransmeridian && b.lng < 0 ? b.lng + M_2PI : b.lng
        sum += (bLng - aLng) * (b.lat + a.lat)
    end
    return sum > 0
end

# These functions depend on higher-layer code and will be defined later:
# _hexRadiusKm, bboxHexEstimate, lineHexEstimate
