"""
    FastH3

Pure Julia implementation of Uber's [H3](https://h3geo.org/) hexagonal hierarchical
geospatial indexing system. Provides the complete H3 v4 public API with no C dependencies.

See the [H3 documentation](https://h3geo.org/docs/) for background on the indexing
system and the [API reference](https://h3geo.org/docs/api/indexing) for the original
C function signatures.
"""
module FastH3

# Layer 0: Types and Constants
include("constants.jl")
include("math_extensions.jl")
include("vec2d.jl")
include("vec3d.jl")
include("lat_lng.jl")

# Layer 1: Core Coordinate Systems
include("coord_ijk.jl")
include("face_ijk.jl")
include("bbox.jl")

# Layer 2: Index and Base Cells
include("base_cells.jl")
include("h3_index.jl")

# Layer 3: Algorithms
include("algos.jl")
include("local_ij.jl")
include("directed_edge.jl")
include("vertex.jl")

# Layer 4: Polygon Operations
include("polygon.jl")
include("linked_geo.jl")
include("area.jl")

# Public API exports
export H3Index, H3Error, H3_NULL, LatLng, CellBoundary, CoordIJ, GeoLoop, GeoPolygon

# Error codes
export E_SUCCESS, E_FAILED, E_DOMAIN, E_LATLNG_DOMAIN, E_RES_DOMAIN,
       E_CELL_INVALID, E_DIR_EDGE_INVALID, E_UNDIR_EDGE_INVALID,
       E_VERTEX_INVALID, E_PENTAGON, E_DUPLICATE_INPUT, E_NOT_NEIGHBORS,
       E_RES_MISMATCH, E_MEMORY_ALLOC, E_MEMORY_BOUNDS, E_OPTION_INVALID

# Containment modes
export ContainmentMode
export CONTAINMENT_CENTER, CONTAINMENT_FULL, CONTAINMENT_OVERLAPPING,
       CONTAINMENT_OVERLAPPING_BBOX, CONTAINMENT_INVALID

# Core functions
export latLngToCell, cellToLatLng, cellToBoundary
export getResolution, getBaseCellNumber, isValidCell, isValidIndex
export isPentagon, isResClassIII
export cellToParent, cellToChildren, cellToChildrenSize
export cellToCenterChild, cellToChildPos, childPosToCell
export compactCells, uncompactCells, uncompactCellsSize
export maxFaceCount, getIcosahedronFaces
export pentagonCount, getPentagons, getRes0Cells

# Grid functions
export gridDisk, gridDiskDistances, gridDiskUnsafe, gridDiskDistancesUnsafe
export gridRingUnsafe, maxGridDiskSize
export gridDistance, gridPathCells, gridPathCellsSize
export cellToLocalIj, localIjToCell
export areNeighborCells

# Directed edge functions
export cellsToDirectedEdge, isValidDirectedEdge
export getDirectedEdgeOrigin, getDirectedEdgeDestination
export directedEdgeToCells, originToDirectedEdges
export directedEdgeToBoundary, reverseDirectedEdge

# Vertex functions
export cellToVertex, cellToVertexes, vertexToLatLng, isValidVertex

# Measurement functions
export greatCircleDistanceRads, greatCircleDistanceKm, greatCircleDistanceM
export degsToRads, radsToDegs
export getHexagonAreaAvgKm2, getHexagonAreaAvgM2
export getHexagonEdgeLengthAvgKm, getHexagonEdgeLengthAvgM
export getNumCells

# String conversion
export stringToH3, h3ToString

# Description
export describeH3Error

end # module FastH3
