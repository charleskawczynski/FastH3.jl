# Geodetic (lat/lng) functions

const EPSILON_DEG = 0.000000001
const EPSILON_RAD = EPSILON_DEG * M_PI_180

@enum LongitudeNormalization::Int32 begin
    NORMALIZE_NONE = 0
    NORMALIZE_EAST = 1
    NORMALIZE_WEST = 2
end

function _posAngleRads(rads::Float64)::Float64
    tmp = rads < 0.0 ? rads + M_2PI : rads
    if rads >= M_2PI
        tmp -= M_2PI
    end
    return tmp
end

function geoAlmostEqualThreshold(p1::LatLng, p2::LatLng, threshold::Float64)::Bool
    return abs(p1.lat - p2.lat) < threshold &&
           abs(p1.lng - p2.lng) < threshold
end

function geoAlmostEqual(p1::LatLng, p2::LatLng)::Bool
    return geoAlmostEqualThreshold(p1, p2, EPSILON_RAD)
end

function setGeoDegs(latDegs::Float64, lngDegs::Float64)::LatLng
    return LatLng(degsToRads(latDegs), degsToRads(lngDegs))
end

degsToRads(degrees::Float64)::Float64 = degrees * M_PI_180
radsToDegs(radians::Float64)::Float64 = radians * M_180_PI

function constrainLat(lat::Float64)::Float64
    while lat > M_PI_2
        lat = lat - M_PI
    end
    return lat
end

function constrainLng(lng::Float64)::Float64
    while lng > M_PI
        lng = lng - 2.0 * M_PI
    end
    while lng < -M_PI
        lng = lng + 2.0 * M_PI
    end
    return lng
end

function normalizeLng(lng::Float64, normalization::LongitudeNormalization)::Float64
    if normalization == NORMALIZE_EAST
        return lng < 0 ? lng + M_2PI : lng
    elseif normalization == NORMALIZE_WEST
        return lng > 0 ? lng - M_2PI : lng
    else
        return lng
    end
end

"""
Great circle distance in radians between two points using the Haversine formula.
"""
function greatCircleDistanceRads(a::LatLng, b::LatLng)::Float64
    sinLat = sin((b.lat - a.lat) * 0.5)
    sinLng = sin((b.lng - a.lng) * 0.5)
    A = sinLat * sinLat + cos(a.lat) * cos(b.lat) * sinLng * sinLng
    return 2.0 * atan(sqrt(A), sqrt(1.0 - A))
end

function greatCircleDistanceKm(a::LatLng, b::LatLng)::Float64
    return greatCircleDistanceRads(a, b) * EARTH_RADIUS_KM
end

function greatCircleDistanceM(a::LatLng, b::LatLng)::Float64
    return greatCircleDistanceKm(a, b) * 1000.0
end

function _geoAzimuthRads(p1::LatLng, p2::LatLng)::Float64
    return atan(cos(p2.lat) * sin(p2.lng - p1.lng),
                cos(p1.lat) * sin(p2.lat) -
                    sin(p1.lat) * cos(p2.lat) * cos(p2.lng - p1.lng))
end

function _geoAzDistanceRads(p1::LatLng, az::Float64, distance::Float64)::LatLng
    if distance < EPSILON
        return p1
    end

    az = _posAngleRads(az)

    if az < EPSILON || abs(az - M_PI) < EPSILON
        if az < EPSILON  # due north
            lat2 = p1.lat + distance
        else  # due south
            lat2 = p1.lat - distance
        end

        if abs(lat2 - M_PI_2) < EPSILON  # north pole
            return LatLng(M_PI_2, 0.0)
        elseif abs(lat2 + M_PI_2) < EPSILON  # south pole
            return LatLng(-M_PI_2, 0.0)
        else
            return LatLng(lat2, constrainLng(p1.lng))
        end
    else
        sinlat = sin(p1.lat) * cos(distance) +
                 cos(p1.lat) * sin(distance) * cos(az)
        sinlat = clamp(sinlat, -1.0, 1.0)
        lat2 = asin(sinlat)

        if abs(lat2 - M_PI_2) < EPSILON  # north pole
            return LatLng(M_PI_2, 0.0)
        elseif abs(lat2 + M_PI_2) < EPSILON  # south pole
            return LatLng(-M_PI_2, 0.0)
        else
            invcosp2lat = 1.0 / cos(lat2)
            sinlng = sin(az) * sin(distance) * invcosp2lat
            coslng = (cos(distance) - sin(p1.lat) * sin(lat2)) /
                     cos(p1.lat) * invcosp2lat
            sinlng = clamp(sinlng, -1.0, 1.0)
            coslng = clamp(coslng, -1.0, 1.0)
            return LatLng(lat2, constrainLng(p1.lng + atan(sinlng, coslng)))
        end
    end
end

# Hexagon area averages (excludes pentagons)
const _hexAreaKm2 = [
    4.357449416078383e+06, 6.097884417941332e+05, 8.680178039899720e+04,
    1.239343465508816e+04, 1.770347654491307e+03, 2.529038581819449e+02,
    3.612906216441245e+01, 5.161293359717191e+00, 7.373275975944177e-01,
    1.053325134272067e-01, 1.504750190766435e-02, 2.149643129451879e-03,
    3.070918756316060e-04, 4.387026794728296e-05, 6.267181135324313e-06,
    8.953115907605790e-07
]

const _hexAreaM2 = [
    4.357449416078390e+12, 6.097884417941339e+11, 8.680178039899731e+10,
    1.239343465508818e+10, 1.770347654491309e+09, 2.529038581819452e+08,
    3.612906216441250e+07, 5.161293359717198e+06, 7.373275975944188e+05,
    1.053325134272069e+05, 1.504750190766437e+04, 2.149643129451882e+03,
    3.070918756316063e+02, 4.387026794728301e+01, 6.267181135324322e+00,
    8.953115907605802e-01
]

const _hexEdgeLengthKm = [
    1281.256011, 483.0568391, 182.5129565, 68.97922179,
    26.07175968, 9.854090990, 3.724532667, 1.406475763,
    0.531414010, 0.200786148, 0.075863783, 0.028663897,
    0.010830188, 0.004092010, 0.001546100, 0.000584169
]

const _hexEdgeLengthM = [
    1281256.011, 483056.8391, 182512.9565, 68979.22179,
    26071.75968, 9854.090990, 3724.532667, 1406.475763,
    531.4140101, 200.7861476, 75.86378287, 28.66389748,
    10.83018784, 4.092010473, 1.546099657, 0.584168630
]

function getHexagonAreaAvgKm2(res::Int)::Tuple{H3Error, Float64}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, 0.0)
    end
    return (E_SUCCESS, _hexAreaKm2[res + 1])
end

function getHexagonAreaAvgM2(res::Int)::Tuple{H3Error, Float64}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, 0.0)
    end
    return (E_SUCCESS, _hexAreaM2[res + 1])
end

function getHexagonEdgeLengthAvgKm(res::Int)::Tuple{H3Error, Float64}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, 0.0)
    end
    return (E_SUCCESS, _hexEdgeLengthKm[res + 1])
end

function getHexagonEdgeLengthAvgM(res::Int)::Tuple{H3Error, Float64}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, 0.0)
    end
    return (E_SUCCESS, _hexEdgeLengthM[res + 1])
end

function getNumCells(res::Int)::Tuple{H3Error, Int64}
    if res < 0 || res > MAX_H3_RES
        return (E_RES_DOMAIN, Int64(0))
    end
    return (E_SUCCESS, Int64(2) + Int64(120) * _ipow(Int64(7), Int64(res)))
end
