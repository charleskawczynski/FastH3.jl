# 3D floating-point vector functions

struct Vec3d
    x::Float64
    y::Float64
    z::Float64
end
Vec3d() = Vec3d(0.0, 0.0, 0.0)

function _pointSquareDist(v1::Vec3d, v2::Vec3d)::Float64
    return (v1.x - v2.x)^2 + (v1.y - v2.y)^2 + (v1.z - v2.z)^2
end

function _geoToVec3d(geo::LatLng)::Vec3d
    r = cos(geo.lat)
    return Vec3d(cos(geo.lng) * r, sin(geo.lng) * r, sin(geo.lat))
end
