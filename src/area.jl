# Area and edge length calculation functions

"""
    cellAreaRads2(h::H3Index) -> (H3Error, Float64)

Compute the exact area of a cell in steradians.

See also the H3 C API: [`cellAreaRads2`](https://h3geo.org/docs/api/misc#cellarearads2)
"""
function cellAreaRads2(h::H3Index)::Tuple{H3Error, Float64}
    err, cb = cellToBoundary(h)
    if err != E_SUCCESS
        return (err, 0.0)
    end
    err2, center = cellToLatLng(h)
    if err2 != E_SUCCESS
        return (err2, 0.0)
    end
    
    area = 0.0
    for i in 1:cb.numVerts
        j = mod1(i + 1, Int(cb.numVerts))
        area += _triangleArea(cb.verts[i], cb.verts[j], center)
    end
    return (E_SUCCESS, area)
end

"""
    cellAreaKm2(h::H3Index) -> (H3Error, Float64)

Compute the exact area of a cell in square kilometers.

See also the H3 C API: [`cellAreaKm2`](https://h3geo.org/docs/api/misc#cellareakm2)
"""
function cellAreaKm2(h::H3Index)::Tuple{H3Error, Float64}
    err, areaRads = cellAreaRads2(h)
    if err != E_SUCCESS
        return (err, 0.0)
    end
    return (E_SUCCESS, areaRads * EARTH_RADIUS_KM * EARTH_RADIUS_KM)
end

"""
    cellAreaM2(h::H3Index) -> (H3Error, Float64)

Compute the exact area of a cell in square meters.

See also the H3 C API: [`cellAreaM2`](https://h3geo.org/docs/api/misc#cellaream2)
"""
function cellAreaM2(h::H3Index)::Tuple{H3Error, Float64}
    err, areaKm = cellAreaKm2(h)
    if err != E_SUCCESS
        return (err, 0.0)
    end
    return (E_SUCCESS, areaKm * 1.0e6)
end

function _triangleArea(a::LatLng, b::LatLng, c::LatLng)::Float64
    # Girard's theorem: area of spherical triangle = |A + B + C - pi|
    az_ab = _geoAzimuthRads(a, b)
    az_ac = _geoAzimuthRads(a, c)
    az_ba = _geoAzimuthRads(b, a)
    az_bc = _geoAzimuthRads(b, c)
    az_ca = _geoAzimuthRads(c, a)
    az_cb = _geoAzimuthRads(c, b)
    
    angle_a = abs(_posAngleRads(az_ab - az_ac))
    if angle_a > M_PI
        angle_a = M_2PI - angle_a
    end
    angle_b = abs(_posAngleRads(az_ba - az_bc))
    if angle_b > M_PI
        angle_b = M_2PI - angle_b
    end
    angle_c = abs(_posAngleRads(az_ca - az_cb))
    if angle_c > M_PI
        angle_c = M_2PI - angle_c
    end
    
    excess = angle_a + angle_b + angle_c - M_PI
    return abs(excess)
end

"""
    edgeLengthRads(edge::H3Index) -> (H3Error, Float64)

Compute the length of a directed edge in radians.

See also the H3 C API: [`edgeLengthRads`](https://h3geo.org/docs/api/misc#edgelengthrads)
"""
function edgeLengthRads(edge::H3Index)::Tuple{H3Error, Float64}
    err, cb = directedEdgeToBoundary(edge)
    if err != E_SUCCESS
        return (err, 0.0)
    end
    length_val = 0.0
    for i in 1:(cb.numVerts - 1)
        length_val += greatCircleDistanceRads(cb.verts[i], cb.verts[i + 1])
    end
    return (E_SUCCESS, length_val)
end

"""
    edgeLengthKm(edge::H3Index) -> (H3Error, Float64)

Compute the length of a directed edge in kilometers.

See also the H3 C API: [`edgeLengthKm`](https://h3geo.org/docs/api/misc#edgelengthkm)
"""
function edgeLengthKm(edge::H3Index)::Tuple{H3Error, Float64}
    err, lengthRads = edgeLengthRads(edge)
    return (err, lengthRads * EARTH_RADIUS_KM)
end

"""
    edgeLengthM(edge::H3Index) -> (H3Error, Float64)

Compute the length of a directed edge in meters.

See also the H3 C API: [`edgeLengthM`](https://h3geo.org/docs/api/misc#edgelengthm)
"""
function edgeLengthM(edge::H3Index)::Tuple{H3Error, Float64}
    err, lengthKm = edgeLengthKm(edge)
    return (err, lengthKm * 1000.0)
end
