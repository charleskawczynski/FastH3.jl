# Constants used throughout the H3 library

const M_PI = 3.14159265358979323846
const M_PI_2 = 1.5707963267948966
const M_2PI = 6.28318530717958647692528676655900576839433
const M_PI_180 = 0.0174532925199432957692369076848861271111
const M_180_PI = 57.29577951308232087679815481410517033240547

const EPSILON = 0.0000000000000001
const M_SQRT3_2 = 0.8660254037844386467637231707529361834714
const M_SIN60 = M_SQRT3_2
const M_RSIN60 = 1.1547005383792515290182975610039149112953

const M_ONETHIRD = 0.333333333333333333333333333333333333333
const M_ONESEVENTH = 0.14285714285714285714285714285714285

const M_AP7_ROT_RADS = 0.333473172251832115336090755351601070065900389
const M_SIN_AP7_ROT = 0.3273268353539885718950318
const M_COS_AP7_ROT = 0.9449111825230680680167902

const EARTH_RADIUS_KM = 6371.007180918475

const RES0_U_GNOMONIC = 0.38196601125010500003
const INV_RES0_U_GNOMONIC = 2.61803398874989588842

const MAX_H3_RES = 15
const NUM_ICOSA_FACES = 20
const NUM_BASE_CELLS = 122
const NUM_HEX_VERTS = 6
const NUM_PENT_VERTS = 5
const NUM_PENTAGONS = 12

const H3_CELL_MODE = 1
const H3_DIRECTEDEDGE_MODE = 2
const H3_EDGE_MODE = 3
const H3_VERTEX_MODE = 4

# H3Index type
"""
64-bit H3 index. Alias for `UInt64`. Encodes cell, directed edge, or vertex identifiers in the H3 system.
"""
const H3Index = UInt64
"""
The null/invalid H3 index (all bits zero).
"""
const H3_NULL = H3Index(0)

# Maximum cell boundary vertices
const MAX_CELL_BNDRY_VERTS = 10

# H3 Error codes
"""
Error codes returned by H3 operations. See also: [`describeH3Error`](@ref).
"""
@enum H3Error::UInt32 begin
    E_SUCCESS = 0
    E_FAILED = 1
    E_DOMAIN = 2
    E_LATLNG_DOMAIN = 3
    E_RES_DOMAIN = 4
    E_CELL_INVALID = 5
    E_DIR_EDGE_INVALID = 6
    E_UNDIR_EDGE_INVALID = 7
    E_VERTEX_INVALID = 8
    E_PENTAGON = 9
    E_DUPLICATE_INPUT = 10
    E_NOT_NEIGHBORS = 11
    E_RES_MISMATCH = 12
    E_MEMORY_ALLOC = 13
    E_MEMORY_BOUNDS = 14
    E_OPTION_INVALID = 15
    E_INDEX_INVALID = 16
    E_BASE_CELL_DOMAIN = 17
    E_DIGIT_DOMAIN = 18
    E_DELETED_DIGIT = 19
end

const H3_ERROR_DESCRIPTIONS = [
    "Success",
    "The operation failed but a more specific error is not available",
    "Argument was outside the acceptable range",
    "Latitude or longitude arguments were outside the acceptable range",
    "Resolution argument was outside the acceptable range",
    "H3Index cell argument was not valid",
    "H3Index directed edge argument was not valid",
    "H3Index undirected edge argument was not valid",
    "H3Index vertex argument was not valid",
    "Pentagon distortion was encountered which the algorithm could not handle",
    "Duplicate input was encountered in the arguments and the algorithm could not handle it",
    "H3Index cell arguments were not neighbors",
    "H3Index cell arguments had incompatible resolutions",
    "Necessary memory allocation failed",
    "Bounds of provided memory were not large enough",
    "Mode or flags argument was not valid",
    "H3Index argument was not valid",
    "Base cell number was outside of acceptable range",
    "Child digits invalid",
    "Deleted subsequence indicates invalid index",
]

"""
    describeH3Error(err::H3Error) -> String

Return a human-readable description of an [`H3Error`](@ref) code.

See also the H3 C API: [`describeH3Error`](https://h3geo.org/docs/api/misc#describeh3error)
"""
function describeH3Error(err::H3Error)
    idx = Int(err) + 1
    if idx < 1 || idx > length(H3_ERROR_DESCRIPTIONS)
        return "Invalid error code"
    end
    return H3_ERROR_DESCRIPTIONS[idx]
end

"""
Mode for polygon containment checks.
"""
@enum ContainmentMode::UInt32 begin
    CONTAINMENT_CENTER = 0
    CONTAINMENT_FULL = 1
    CONTAINMENT_OVERLAPPING = 2
    CONTAINMENT_OVERLAPPING_BBOX = 3
    CONTAINMENT_INVALID = 4
end

# Core geometric types

"""
Latitude/longitude coordinate pair in radians.
"""
struct LatLng
    lat::Float64
    lng::Float64
end
LatLng() = LatLng(0.0, 0.0)

"""
Cell boundary represented as a fixed-capacity array of [`LatLng`](@ref) vertices.
"""
struct CellBoundary
    numVerts::Int32
    verts::NTuple{MAX_CELL_BNDRY_VERTS, LatLng}
end
function CellBoundary()
    CellBoundary(Int32(0), ntuple(_ -> LatLng(), MAX_CELL_BNDRY_VERTS))
end

"""
Ordered sequence of [`LatLng`](@ref) vertices forming a closed loop.
"""
struct GeoLoop
    numVerts::Int32
    verts::Vector{LatLng}
end
GeoLoop() = GeoLoop(Int32(0), LatLng[])

"""
Polygon defined by an outer [`GeoLoop`](@ref) and zero or more holes.
"""
struct GeoPolygon
    geoloop::GeoLoop
    numHoles::Int32
    holes::Vector{GeoLoop}
end
GeoPolygon() = GeoPolygon(GeoLoop(), Int32(0), GeoLoop[])

struct GeoMultiPolygon
    numPolygons::Int32
    polygons::Vector{GeoPolygon}
end
GeoMultiPolygon() = GeoMultiPolygon(Int32(0), GeoPolygon[])

"""
IJ coordinate pair used for local grid coordinates. See the H3 C API: [`CoordIJ`](https://h3geo.org/docs/api/traversal).
"""
struct CoordIJ
    i::Int32
    j::Int32
end
CoordIJ() = CoordIJ(Int32(0), Int32(0))

# Linked geo structures (Vector-based replacements for C linked lists)
struct LinkedGeoLoop
    verts::Vector{LatLng}
end
LinkedGeoLoop() = LinkedGeoLoop(LatLng[])

struct LinkedGeoPolygon
    loops::Vector{LinkedGeoLoop}
end
LinkedGeoPolygon() = LinkedGeoPolygon(LinkedGeoLoop[])
