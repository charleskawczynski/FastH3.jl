# 2D floating-point vector functions

struct Vec2d
    x::Float64
    y::Float64
end
Vec2d() = Vec2d(0.0, 0.0)

function _v2dMag(v::Vec2d)::Float64
    return sqrt(v.x * v.x + v.y * v.y)
end

"""
Find the intersection between two lines. Assumes that the lines intersect
and that the intersection is not at an endpoint of either line.
"""
function _v2dIntersect(p0::Vec2d, p1::Vec2d, p2::Vec2d, p3::Vec2d)::Vec2d
    s1x = p1.x - p0.x
    s1y = p1.y - p0.y
    s2x = p3.x - p2.x
    s2y = p3.y - p2.y

    t = (s2x * (p0.y - p2.y) - s2y * (p0.x - p2.x)) /
        (-s2x * s1y + s1x * s2y)

    return Vec2d(p0.x + t * s1x, p0.y + t * s1y)
end

function _v2dAlmostEquals(v1::Vec2d, v2::Vec2d)::Bool
    return abs(v1.x - v2.x) < eps(Float32) &&
           abs(v1.y - v2.y) < eps(Float32)
end
