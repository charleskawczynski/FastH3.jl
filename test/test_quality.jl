using Aqua
using JET
using CodeComplexity

function _qc_noop_h3(::H3Index) end
function _qc_noop_h3_ring(::H3Index, ::Int) end
function _qc_noop_int(::Int) end

@testset "Code quality" begin

    @testset "Aqua" begin
        Aqua.test_all(FastH3)
    end

    @testset "JET optimization" begin
        cell = H3Index(0x088283080ddfffff)
        cell2 = H3Index(0x088283080dc3ffff)
        latlng = LatLng(0.0, 0.0)
        latlng2 = LatLng(0.1, 0.1)
        ij = CoordIJ(Int32(0), Int32(0))
        cells = H3Index[cell]

        @testset "Index inspection" begin
            @test_opt getResolution(cell)
            @test_opt getBaseCellNumber(cell)
            @test_opt isValidCell(cell)
            @test_opt isPentagon(cell)
            @test_opt isResClassIII(cell)
            @test_opt isValidIndex(cell)
        end

        @testset "Coordinate conversion" begin
            @test_opt latLngToCell(latlng, 5)
            @test_opt cellToLatLng(cell)
            @test_opt cellToBoundary(cell)
        end

        @testset "Hierarchy" begin
            @test_opt cellToParent(cell, 5)
            @test_opt cellToCenterChild(cell, 10)
            @test_opt cellToChildrenSize(cell, 10)
            @test_opt cellToChildren(cell, 10)
            @test_opt cellToChildren(_qc_noop_h3, cell, 10)
        end

        @testset "Compact / uncompact" begin
            @test_opt compactCells(cells)
            @test_opt uncompactCellsSize(cells, 10)
            @test_opt uncompactCells(cells, 10)
            @test_opt uncompactCells(_qc_noop_h3, cells, 10)
        end

        @testset "Grid traversal" begin
            @test_opt maxGridDiskSize(3)
            @test_opt gridDisk(cell, 1)
            @test_opt gridDiskDistances(cell, 1)
            @test_opt gridDiskUnsafe(cell, 1)
            @test_opt gridDiskDistancesUnsafe(cell, 1)
            @test_opt gridRingUnsafe(cell, 1)
            @test_opt gridDiskUnsafe(_qc_noop_h3, cell, 1)
            @test_opt gridDiskDistancesUnsafe(_qc_noop_h3_ring, cell, 1)
            @test_opt gridRingUnsafe(_qc_noop_h3, cell, 1)
            @test_opt gridDistance(cell, cell2)
            @test_opt gridPathCellsSize(cell, cell2)
            @test_opt gridPathCells(cell, cell2)
            @test_opt gridPathCells(_qc_noop_h3, cell, cell2)
            @test_opt cellToLocalIj(cell, cell2, UInt32(0))
            @test_opt localIjToCell(cell, ij, UInt32(0))
            @test_opt areNeighborCells(cell, cell2)
        end

        @testset "Directed edges" begin
            @test_opt isValidDirectedEdge(cell)
            @test_opt getDirectedEdgeOrigin(cell)
            @test_opt getDirectedEdgeDestination(cell)
            @test_opt directedEdgeToCells(cell)
            @test_opt cellsToDirectedEdge(cell, cell2)
            @test_opt originToDirectedEdges(cell)
            @test_opt originToDirectedEdges(_qc_noop_h3, cell)
            @test_opt directedEdgeToBoundary(cell)
            @test_opt reverseDirectedEdge(cell)
        end

        @testset "Vertexes" begin
            @test_opt cellToVertex(cell, 0)
            @test_opt cellToVertexes(cell)
            @test_opt cellToVertexes(_qc_noop_h3, cell)
            @test_opt vertexToLatLng(cell)
            @test_opt isValidVertex(cell)
        end

        @testset "Measurement" begin
            @test_opt greatCircleDistanceRads(latlng, latlng2)
            @test_opt greatCircleDistanceKm(latlng, latlng2)
            @test_opt greatCircleDistanceM(latlng, latlng2)
            @test_opt degsToRads(1.0)
            @test_opt radsToDegs(1.0)
            @test_opt getHexagonAreaAvgKm2(5)
            @test_opt getHexagonAreaAvgM2(5)
            @test_opt getHexagonEdgeLengthAvgKm(5)
            @test_opt getHexagonEdgeLengthAvgM(5)
            @test_opt getNumCells(5)
            @test_opt pentagonCount()
        end

        @testset "Faces" begin
            @test_opt maxFaceCount(cell)
            @test_opt getIcosahedronFaces(cell)
            @test_opt getIcosahedronFaces(_qc_noop_int, cell)
        end

        @testset "Pentagons and res0" begin
            @test_opt getPentagons(5)
            @test_opt getPentagons(_qc_noop_h3, 5)
            @test_opt getRes0Cells()
            @test_opt getRes0Cells(_qc_noop_h3)
        end

        @testset "String conversion" begin
            if VERSION ≥ v"1.11"
                @test_opt stringToH3("088283080ddfffff")
            else
                @test_opt broken=true stringToH3("088283080ddfffff")
            end
            @test_opt h3ToString(cell)
            @test_opt describeH3Error(E_SUCCESS)
        end
    end

    @testset "CodeComplexity" begin
        check_complexity(FastH3; max_complexity=40)
        @test true
    end

    @testset "Allocations" begin
        latlng = LatLng(0.0, 0.0)
        latlng2 = LatLng(0.1, 0.1)
        err, cell = latLngToCell(latlng, 5)
        @assert err == E_SUCCESS
        err, disk = gridDisk(cell, 1)
        @assert err == E_SUCCESS
        neighbor = disk[2]
        err, edge = cellsToDirectedEdge(cell, neighbor)
        @assert err == E_SUCCESS
        err, ij = cellToLocalIj(cell, neighbor, UInt32(0))
        @assert err == E_SUCCESS
        cells_cb = H3Index[cell]

        # Warmup
        getResolution(cell)
        getBaseCellNumber(cell)
        isValidCell(cell)
        isPentagon(cell)
        isResClassIII(cell)
        isValidDirectedEdge(edge)
        isValidVertex(cell)
        degsToRads(1.0)
        radsToDegs(1.0)
        greatCircleDistanceRads(latlng, latlng2)
        greatCircleDistanceKm(latlng, latlng2)
        greatCircleDistanceM(latlng, latlng2)
        maxGridDiskSize(3)
        getHexagonAreaAvgKm2(5)
        getHexagonAreaAvgM2(5)
        getHexagonEdgeLengthAvgKm(5)
        getHexagonEdgeLengthAvgM(5)
        getNumCells(5)
        pentagonCount()
        describeH3Error(E_SUCCESS)
        latLngToCell(latlng, 5)
        cellToLatLng(cell)
        cellToParent(cell, 3)
        cellToCenterChild(cell, 8)
        cellToChildrenSize(cell, 8)
        getDirectedEdgeOrigin(edge)
        getDirectedEdgeDestination(edge)
        reverseDirectedEdge(edge)
        directedEdgeToCells(edge)
        areNeighborCells(cell, neighbor)
        gridDistance(cell, neighbor)
        cellToLocalIj(cell, neighbor, UInt32(0))
        localIjToCell(cell, ij, UInt32(0))
        gridPathCellsSize(cell, neighbor)
        cellsToDirectedEdge(cell, neighbor)
        cellToChildren(_qc_noop_h3, cell, 8)
        uncompactCells(_qc_noop_h3, cells_cb, 8)
        getPentagons(_qc_noop_h3, 5)
        getRes0Cells(_qc_noop_h3)
        getIcosahedronFaces(_qc_noop_int, cell)
        gridDiskDistancesUnsafe(_qc_noop_h3_ring, cell, 1)
        gridDiskUnsafe(_qc_noop_h3, cell, 1)
        gridRingUnsafe(_qc_noop_h3, cell, 1)
        gridPathCells(_qc_noop_h3, cell, neighbor)
        cellToVertexes(_qc_noop_h3, cell)
        originToDirectedEdges(_qc_noop_h3, cell)

        @testset "Bit manipulation" begin
            @test @allocated(getResolution(cell)) == 0
            @test @allocated(getBaseCellNumber(cell)) == 0
            @test @allocated(isValidCell(cell)) == 0
            @test @allocated(isPentagon(cell)) == 0
            @test @allocated(isResClassIII(cell)) == 0
            @test @allocated(isValidDirectedEdge(edge)) == 0
            @test @allocated(isValidVertex(cell)) == 0
        end

        @testset "Distance and math" begin
            @test @allocated(degsToRads(1.0)) == 0
            @test @allocated(radsToDegs(1.0)) == 0
            @test @allocated(greatCircleDistanceRads(latlng, latlng2)) == 0
            @test @allocated(greatCircleDistanceKm(latlng, latlng2)) == 0
            @test @allocated(greatCircleDistanceM(latlng, latlng2)) == 0
        end

        @testset "Table lookups" begin
            @test @allocated(maxGridDiskSize(3)) == 0
            @test @allocated(getHexagonAreaAvgKm2(5)) == 0
            @test @allocated(getHexagonAreaAvgM2(5)) == 0
            @test @allocated(getHexagonEdgeLengthAvgKm(5)) == 0
            @test @allocated(getHexagonEdgeLengthAvgM(5)) == 0
            @test @allocated(getNumCells(5)) == 0
            @test @allocated(pentagonCount()) == 0
        end

        @testset "Core index operations" begin
            @test @allocated(latLngToCell(latlng, 5)) == 0
            @test @allocated(cellToLatLng(cell)) == 0
            @test @allocated(cellToParent(cell, 3)) == 0
            @test @allocated(cellToCenterChild(cell, 8)) == 0
            @test @allocated(cellToChildrenSize(cell, 8)) == 0
        end

        @testset "Directed edge operations" begin
            @test @allocated(getDirectedEdgeOrigin(edge)) == 0
            @test @allocated(getDirectedEdgeDestination(edge)) == 0
            @test @allocated(reverseDirectedEdge(edge)) == 0
            @test @allocated(directedEdgeToCells(edge)) == 0
        end

        @testset "Grid scalar operations" begin
            @test @allocated(areNeighborCells(cell, neighbor)) == 0
            @test @allocated(gridDistance(cell, neighbor)) == 0
            @test @allocated(cellToLocalIj(cell, neighbor, UInt32(0))) == 0
            @test @allocated(localIjToCell(cell, ij, UInt32(0))) == 0
            @test @allocated(gridPathCellsSize(cell, neighbor)) == 0
            @test @allocated(cellsToDirectedEdge(cell, neighbor)) == 0
        end

        @testset "Callback enumeration (no heap for output vectors)" begin
            @test @allocated(cellToChildren(_qc_noop_h3, cell, 8)) == 0
            @test @allocated(uncompactCells(_qc_noop_h3, cells_cb, 8)) == 0
            @test @allocated(getPentagons(_qc_noop_h3, 5)) == 0
            @test @allocated(getRes0Cells(_qc_noop_h3)) == 0
            @test @allocated(getIcosahedronFaces(_qc_noop_int, cell)) == 0
            @test @allocated(gridDiskDistancesUnsafe(_qc_noop_h3_ring, cell, 1)) == 0
            @test @allocated(gridDiskUnsafe(_qc_noop_h3, cell, 1)) == 0
            @test @allocated(gridRingUnsafe(_qc_noop_h3, cell, 1)) == 0
            @test @allocated(gridPathCells(_qc_noop_h3, cell, neighbor)) == 0
            @test @allocated(cellToVertexes(_qc_noop_h3, cell)) == 0
            @test @allocated(originToDirectedEdges(_qc_noop_h3, cell)) == 0
        end
    end
end
