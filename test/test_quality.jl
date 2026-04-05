using Aqua
using JET
using CodeComplexity

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
            @test_opt target_modules = (FastH3,) getResolution(cell)
            @test_opt target_modules = (FastH3,) getBaseCellNumber(cell)
            @test_opt target_modules = (FastH3,) isValidCell(cell)
            @test_opt target_modules = (FastH3,) isPentagon(cell)
            @test_opt target_modules = (FastH3,) isResClassIII(cell)
            @test_opt target_modules = (FastH3,) isValidIndex(cell)
        end

        @testset "Coordinate conversion" begin
            @test_opt target_modules = (FastH3,) latLngToCell(latlng, 5)
            @test_opt target_modules = (FastH3,) cellToLatLng(cell)
            @test_opt target_modules = (FastH3,) cellToBoundary(cell)
        end

        @testset "Hierarchy" begin
            @test_opt target_modules = (FastH3,) cellToParent(cell, 5)
            @test_opt target_modules = (FastH3,) cellToCenterChild(cell, 10)
            @test_opt target_modules = (FastH3,) cellToChildrenSize(cell, 10)
            @test_opt target_modules = (FastH3,) cellToChildren(cell, 10)
        end

        @testset "Compact / uncompact" begin
            @test_opt target_modules = (FastH3,) compactCells(cells)
            @test_opt target_modules = (FastH3,) uncompactCellsSize(cells, 10)
            @test_opt target_modules = (FastH3,) uncompactCells(cells, 10)
        end

        @testset "Grid traversal" begin
            @test_opt target_modules = (FastH3,) maxGridDiskSize(3)
            @test_opt target_modules = (FastH3,) gridDisk(cell, 1)
            @test_opt target_modules = (FastH3,) gridDiskDistances(cell, 1)
            @test_opt target_modules = (FastH3,) gridDiskUnsafe(cell, 1)
            @test_opt target_modules = (FastH3,) gridDiskDistancesUnsafe(cell, 1)
            @test_opt target_modules = (FastH3,) gridRingUnsafe(cell, 1)
            @test_opt target_modules = (FastH3,) gridDistance(cell, cell2)
            @test_opt target_modules = (FastH3,) gridPathCellsSize(cell, cell2)
            @test_opt target_modules = (FastH3,) gridPathCells(cell, cell2)
            @test_opt target_modules = (FastH3,) cellToLocalIj(cell, cell2, UInt32(0))
            @test_opt target_modules = (FastH3,) localIjToCell(cell, ij, UInt32(0))
            @test_opt target_modules = (FastH3,) areNeighborCells(cell, cell2)
        end

        @testset "Directed edges" begin
            @test_opt target_modules = (FastH3,) isValidDirectedEdge(cell)
            @test_opt target_modules = (FastH3,) getDirectedEdgeOrigin(cell)
            @test_opt target_modules = (FastH3,) getDirectedEdgeDestination(cell)
            @test_opt target_modules = (FastH3,) directedEdgeToCells(cell)
            @test_opt target_modules = (FastH3,) cellsToDirectedEdge(cell, cell2)
            @test_opt target_modules = (FastH3,) originToDirectedEdges(cell)
            @test_opt target_modules = (FastH3,) directedEdgeToBoundary(cell)
            @test_opt target_modules = (FastH3,) reverseDirectedEdge(cell)
        end

        @testset "Vertexes" begin
            @test_opt target_modules = (FastH3,) cellToVertex(cell, 0)
            @test_opt target_modules = (FastH3,) cellToVertexes(cell)
            @test_opt target_modules = (FastH3,) vertexToLatLng(cell)
            @test_opt target_modules = (FastH3,) isValidVertex(cell)
        end

        @testset "Measurement" begin
            @test_opt target_modules = (FastH3,) greatCircleDistanceRads(latlng, latlng2)
            @test_opt target_modules = (FastH3,) greatCircleDistanceKm(latlng, latlng2)
            @test_opt target_modules = (FastH3,) greatCircleDistanceM(latlng, latlng2)
            @test_opt target_modules = (FastH3,) degsToRads(1.0)
            @test_opt target_modules = (FastH3,) radsToDegs(1.0)
            @test_opt target_modules = (FastH3,) getHexagonAreaAvgKm2(5)
            @test_opt target_modules = (FastH3,) getHexagonAreaAvgM2(5)
            @test_opt target_modules = (FastH3,) getHexagonEdgeLengthAvgKm(5)
            @test_opt target_modules = (FastH3,) getHexagonEdgeLengthAvgM(5)
            @test_opt target_modules = (FastH3,) getNumCells(5)
            @test_opt target_modules = (FastH3,) pentagonCount()
        end

        @testset "Faces" begin
            @test_opt target_modules = (FastH3,) maxFaceCount(cell)
            @test_opt target_modules = (FastH3,) getIcosahedronFaces(cell)
        end

        @testset "Pentagons and res0" begin
            @test_opt target_modules = (FastH3,) getPentagons(5)
            @test_opt target_modules = (FastH3,) getRes0Cells()
        end

        @testset "String conversion" begin
            @test_opt target_modules = (FastH3,) stringToH3("088283080ddfffff")
            @test_opt target_modules = (FastH3,) h3ToString(cell)
            @test_opt target_modules = (FastH3,) describeH3Error(E_SUCCESS)
        end
    end

    @testset "CodeComplexity" begin
        check_complexity(FastH3; max_complexity = 40)
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
    end
end
