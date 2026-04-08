module FastH3Extension

import ..FastH3: H3Index, H3_NULL, H3Error, LatLng, CellBoundary,
    E_SUCCESS, E_RES_MISMATCH, E_CELL_INVALID,
    latLngToCell, cellToLatLng, cellToBoundary,
    areNeighborCells, getResolution,
    greatCircleDistanceRads, EARTH_RADIUS_KM,
    getHexagonEdgeLengthAvgKm,
    getBaseCellNumber,
    gridDisk, gridDistance, isValidCell

import ..FastH3: gridPathCells as _cube_gridPathCells

include("great_circle_path.jl")

end # module FastH3Extension
