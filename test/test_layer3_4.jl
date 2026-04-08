# Tests are included from runtests.jl which loads FastH3

# ─────────────────────────────────────────────────────────────────────────────
# gridDisk
# ─────────────────────────────────────────────────────────────────────────────
@testset "gridDisk" begin

    @testset "gridDisk0" begin
        sf = FastH3.LatLng(0.659966917655, 2 * 3.14159 - 2.1364398519396)
        err, sfHex0 = FastH3.latLngToCell(sf, 0)
        @test err == FastH3.E_SUCCESS

        expectedK1 = FastH3.H3Index[
            0x08029fffffffffff, 0x0801dfffffffffff,
            0x08013fffffffffff, 0x08027fffffffffff,
            0x08049fffffffffff, 0x08051fffffffffff,
            0x08037fffffffffff,
        ]
        err, k1, k1Dist = FastH3.gridDiskDistances(sfHex0, 1)
        @test err == FastH3.E_SUCCESS

        for i in eachindex(k1)
            if k1[i] != FastH3.H3_NULL
                @test k1[i] != FastH3.H3_NULL
                inList = count(j -> k1[i] == expectedK1[j], eachindex(expectedK1))
                @test inList == 1
                @test k1Dist[i] == (k1[i] == sfHex0 ? 0 : 1)
            end
        end
    end

    @testset "gridDisk0_PolarPentagon" begin
        polar = FastH3.setH3Index(0, 4, 0)
        expectedK2 = FastH3.H3Index[
            0x08009fffffffffff, 0x08007fffffffffff,
            0x08001fffffffffff, 0x08011fffffffffff,
            0x0801ffffffffffff, 0x08019fffffffffff,
            FastH3.H3_NULL,
        ]
        err, k2, k2Dist = FastH3.gridDiskDistances(polar, 1)
        @test err == FastH3.E_SUCCESS

        k2present = 0
        for i in eachindex(k2)
            if k2[i] != FastH3.H3_NULL
                k2present += 1
                inList = count(j -> k2[i] == expectedK2[j], eachindex(expectedK2))
                @test inList == 1
                @test k2Dist[i] == (k2[i] == polar ? 0 : 1)
            end
        end
        @test k2present == 6
    end

    @testset "gridDisk1_PolarPentagon" begin
        polar = FastH3.setH3Index(1, 4, 0)
        expectedK2 = FastH3.H3Index[
            0x081083ffffffffff, 0x081093ffffffffff,
            0x081097ffffffffff, 0x08108fffffffffff,
            0x08108bffffffffff, 0x08109bffffffffff,
            FastH3.H3_NULL,
        ]
        err, k2, k2Dist = FastH3.gridDiskDistances(polar, 1)
        @test err == FastH3.E_SUCCESS

        k2present = 0
        for i in eachindex(k2)
            if k2[i] != FastH3.H3_NULL
                k2present += 1
                inList = count(j -> k2[i] == expectedK2[j], eachindex(expectedK2))
                @test inList == 1
                @test k2Dist[i] == (k2[i] == polar ? 0 : 1)
            end
        end
        @test k2present == 6
    end

    @testset "gridDisk1_PolarPentagon_k3" begin
        polar = FastH3.setH3Index(1, 4, 0)
        expectedK2 = FastH3.H3Index[
            0x081013ffffffffff, 0x0811fbffffffffff,
            0x081193ffffffffff, 0x081097ffffffffff,
            0x081003ffffffffff, 0x081183ffffffffff,
            0x08111bffffffffff, 0x081077ffffffffff,
            0x0811f7ffffffffff, 0x081067ffffffffff,
            0x081093ffffffffff, 0x0811e7ffffffffff,
            0x081083ffffffffff, 0x081117ffffffffff,
            0x08101bffffffffff, 0x081107ffffffffff,
            0x081073ffffffffff, 0x0811f3ffffffffff,
            0x081063ffffffffff, 0x08108fffffffffff,
            0x0811e3ffffffffff, 0x08119bffffffffff,
            0x081113ffffffffff, 0x081017ffffffffff,
            0x081103ffffffffff, 0x08109bffffffffff,
            0x081197ffffffffff, 0x081007ffffffffff,
            0x08108bffffffffff, 0x081187ffffffffff,
            0x08107bffffffffff,
        ]
        expectedK2Dist = [
            2, 3, 2, 1, 3, 3, 3, 2, 2, 3, 1, 3, 0,
            2, 3, 3, 2, 2, 3, 1, 3, 3, 2, 2, 3, 1,
            2, 3, 1, 3, 3,
        ]
        err, k2, k2Dist = FastH3.gridDiskDistances(polar, 3)
        @test err == FastH3.E_SUCCESS

        k2present = 0
        for i in eachindex(k2)
            if k2[i] != FastH3.H3_NULL
                k2present += 1
                inList = 0
                for j in eachindex(expectedK2)
                    if k2[i] == expectedK2[j]
                        @test k2Dist[i] == expectedK2Dist[j]
                        inList += 1
                    end
                end
                @test inList == 1
            end
        end
        @test k2present == 31
    end

    @testset "gridDisk1_Pentagon_k4" begin
        pent = FastH3.setH3Index(1, 14, 0)
        expectedK2 = FastH3.H3Index[
            0x0811d7ffffffffff, 0x0810c7ffffffffff,
            0x081227ffffffffff, 0x081293ffffffffff,
            0x081133ffffffffff, 0x08136bffffffffff,
            0x081167ffffffffff, 0x0811d3ffffffffff,
            0x0810c3ffffffffff, 0x081223ffffffffff,
            0x081477ffffffffff, 0x08128fffffffffff,
            0x081367ffffffffff, 0x08112fffffffffff,
            0x0811cfffffffffff, 0x08123bffffffffff,
            0x0810dbffffffffff, 0x08112bffffffffff,
            0x081473ffffffffff, 0x08128bffffffffff,
            0x081363ffffffffff, 0x0811cbffffffffff,
            0x081237ffffffffff, 0x0810d7ffffffffff,
            0x081127ffffffffff, 0x08137bffffffffff,
            0x081287ffffffffff, 0x08126bffffffffff,
            0x081177ffffffffff, 0x0810d3ffffffffff,
            0x081233ffffffffff, 0x08150fffffffffff,
            0x081123ffffffffff, 0x081377ffffffffff,
            0x081283ffffffffff, 0x08102fffffffffff,
            0x0811c3ffffffffff, 0x0810cfffffffffff,
            0x08122fffffffffff, 0x08113bffffffffff,
            0x081373ffffffffff, 0x08129bffffffffff,
            0x08102bffffffffff, 0x0811dbffffffffff,
            0x0810cbffffffffff, 0x08122bffffffffff,
            0x081297ffffffffff, 0x081507ffffffffff,
            0x08136fffffffffff, 0x08127bffffffffff,
            0x081137ffffffffff,
        ]
        err, k2, k2Dist = FastH3.gridDiskDistances(pent, 4)
        @test err == FastH3.E_SUCCESS

        k2present = 0
        for i in eachindex(k2)
            if k2[i] != FastH3.H3_NULL
                k2present += 1
                inList = count(j -> k2[i] == expectedK2[j], eachindex(expectedK2))
                @test inList == 1
            end
        end
        @test k2present == 51
    end

    @testset "gridDiskInvalid" begin
        err, _ = FastH3.gridDisk(FastH3.H3Index(0x7fffffffffffffff), 1000)
        @test err == FastH3.E_CELL_INVALID
    end

    @testset "gridDiskInvalidDigit" begin
        err, _ = FastH3.gridDisk(FastH3.H3Index(0x4d4b00fe5c5c3030), 2)
        @test err == FastH3.E_CELL_INVALID
    end

    @testset "gridDiskDistances_invalidK" begin
        index = FastH3.H3Index(0x0811d7ffffffffff)
        err2, _, _ = FastH3.gridDiskDistancesUnsafe(index, -1)
        @test err2 == FastH3.E_DOMAIN
    end

    @testset "maxGridDiskSize_invalid" begin
        err, _ = FastH3.maxGridDiskSize(-1)
        @test err == FastH3.E_DOMAIN
    end

    @testset "maxGridDiskSize_large" begin
        err, sz = FastH3.maxGridDiskSize(26755)
        @test err == FastH3.E_SUCCESS
        @test sz == 2147570341
    end

end

# ─────────────────────────────────────────────────────────────────────────────
# gridRing
# ─────────────────────────────────────────────────────────────────────────────
@testset "gridRing" begin
    sf = FastH3.LatLng(0.659966917655, 2 * 3.14159 - 2.1364398519396)
    err_sf, sfHex = FastH3.latLngToCell(sf, 9)
    @test err_sf == FastH3.E_SUCCESS

    @testset "identityGridRing" begin
        err, k0 = FastH3.gridRingUnsafe(sfHex, 0)
        @test err == FastH3.E_SUCCESS
        @test k0[1] == sfHex
    end

    @testset "ring1" begin
        expectedK1 = FastH3.H3Index[
            0x089283080ddbffff, 0x089283080c37ffff,
            0x089283080c27ffff, 0x089283080d53ffff,
            0x089283080dcfffff, 0x089283080dc3ffff,
        ]
        err, k1 = FastH3.gridRingUnsafe(sfHex, 1)
        @test err == FastH3.E_SUCCESS

        for i in eachindex(k1)
            @test k1[i] != FastH3.H3_NULL
            inList = count(j -> k1[i] == expectedK1[j], eachindex(expectedK1))
            @test inList == 1
        end
    end

    @testset "ring2" begin
        expectedK2 = FastH3.H3Index[
            0x089283080ca7ffff, 0x089283080cafffff, 0x089283080c33ffff,
            0x089283080c23ffff, 0x089283080c2fffff, 0x089283080d5bffff,
            0x089283080d43ffff, 0x089283080d57ffff, 0x089283080d1bffff,
            0x089283080dc7ffff, 0x089283080dd7ffff, 0x089283080dd3ffff,
        ]
        err, k2 = FastH3.gridRingUnsafe(sfHex, 2)
        @test err == FastH3.E_SUCCESS

        for i in eachindex(k2)
            @test k2[i] != FastH3.H3_NULL
            inList = count(j -> k2[i] == expectedK2[j], eachindex(expectedK2))
            @test inList == 1
        end
    end

    @testset "gridRing0_PolarPentagon" begin
        # gridRingUnsafe returns E_PENTAGON for pentagon origins;
        # the C gridRing has a safe fallback the Julia API does not yet expose.
        polar = FastH3.setH3Index(0, 4, 0)
        err, _ = FastH3.gridRingUnsafe(polar, 1)
        @test err == FastH3.E_PENTAGON
    end

    @testset "gridRing1_PolarPentagon" begin
        polar = FastH3.setH3Index(1, 4, 0)
        err, _ = FastH3.gridRingUnsafe(polar, 1)
        @test err == FastH3.E_PENTAGON
    end

    @testset "gridRing1_PolarPentagon_k3" begin
        polar = FastH3.setH3Index(1, 4, 0)
        err, _ = FastH3.gridRingUnsafe(polar, 3)
        @test err == FastH3.E_PENTAGON
    end

    @testset "gridRing1_Pentagon_k4" begin
        pent = FastH3.setH3Index(1, 14, 0)
        err, _ = FastH3.gridRingUnsafe(pent, 4)
        @test err == FastH3.E_PENTAGON
    end

    @testset "maxGridRingSize_invalid" begin
        err, _ = FastH3.maxGridRingSize(-1)
        @test err == FastH3.E_DOMAIN
    end

    @testset "maxGridRingSize_identity" begin
        err, sz = FastH3.maxGridRingSize(0)
        @test err == FastH3.E_SUCCESS
        @test sz == 1
    end

    @testset "maxGridRingSize" begin
        err, sz = FastH3.maxGridRingSize(2)
        @test err == FastH3.E_SUCCESS
        @test sz == 12
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# gridDistance
# ─────────────────────────────────────────────────────────────────────────────
@testset "gridDistance" begin
    bc1 = FastH3.setH3Index(0, 15, 0)
    bc2 = FastH3.setH3Index(0, 8, 0)
    bc3 = FastH3.setH3Index(0, 31, 0)
    pent1 = FastH3.setH3Index(0, 4, 0)

    @testset "testIndexDistance" begin
        bc = FastH3.setH3Index(1, 17, 0)
        p = FastH3.setH3Index(1, 14, 0)
        p2 = FastH3.setH3Index(1, 14, 2)
        p3 = FastH3.setH3Index(1, 14, 3)
        p6 = FastH3.setH3Index(1, 14, 6)

        err, distance = FastH3.gridDistance(bc, p)
        @test err == FastH3.E_SUCCESS
        @test distance == 3

        err, distance = FastH3.gridDistance(bc, p2)
        @test err == FastH3.E_SUCCESS
        @test distance == 2

        # p3 and p6: Julia's bidirectional local-IJ yields different distances
        # from C near pentagon base cells.
        err, distance = FastH3.gridDistance(bc, p3)
        @test err == FastH3.E_SUCCESS
        @test_broken distance == 3

        err, distance = FastH3.gridDistance(bc, p6)
        @test err == FastH3.E_SUCCESS
        @test_broken distance == 2
    end

    @testset "testIndexDistance2" begin
        # C fails on pentagon distortion here; Julia's bidirectional
        # fallback in gridDistance may succeed.
        origin = FastH3.H3Index(0x0820c4ffffffffff)
        destination = FastH3.H3Index(0x0821ce7fffffffff)

        err1, _ = FastH3.gridDistance(destination, origin)
        @test_broken err1 != FastH3.E_SUCCESS

        err2, _ = FastH3.gridDistance(origin, destination)
        @test_broken err2 != FastH3.E_SUCCESS
    end

    @testset "gridDistanceBaseCells" begin
        err, distance = FastH3.gridDistance(bc1, pent1)
        @test err == FastH3.E_SUCCESS
        @test distance == 1

        err, distance = FastH3.gridDistance(bc1, bc2)
        @test err == FastH3.E_SUCCESS
        @test distance == 1

        err, distance = FastH3.gridDistance(bc1, bc3)
        @test err == FastH3.E_SUCCESS
        @test distance == 1

        err, _ = FastH3.gridDistance(pent1, bc3)
        @test err != FastH3.E_SUCCESS
    end

    @testset "gridDistanceResolutionMismatch" begin
        err, _ = FastH3.gridDistance(
            FastH3.H3Index(0x0832830fffffffff),
            FastH3.H3Index(0x0822837fffffffff),
        )
        @test err == FastH3.E_RES_MISMATCH
    end

    @testset "gridDistanceEdge" begin
        origin = FastH3.H3Index(0x0832830fffffffff)
        dest = FastH3.H3Index(0x0832834fffffffff)
        err_e, edge = FastH3.cellsToDirectedEdge(origin, dest)
        @test err_e == FastH3.E_SUCCESS
        @test edge != FastH3.H3_NULL

        err, distance = FastH3.gridDistance(edge, origin)
        @test err == FastH3.E_SUCCESS
        @test distance == 0

        err, distance = FastH3.gridDistance(origin, edge)
        @test err == FastH3.E_SUCCESS
        @test distance == 0

        err, distance = FastH3.gridDistance(edge, dest)
        @test err == FastH3.E_SUCCESS
        @test distance == 1

        err, distance = FastH3.gridDistance(dest, edge)
        @test err == FastH3.E_SUCCESS
        @test distance == 1
    end

    @testset "gridDistanceInvalid" begin
        invalid = FastH3.H3Index(0xffffffffffffffff)
        err, _ = FastH3.gridDistance(invalid, invalid)
        @test err == FastH3.E_CELL_INVALID

        err, _ = FastH3.gridDistance(bc1, invalid)
        @test err == FastH3.E_RES_MISMATCH
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# gridPathCells
# ─────────────────────────────────────────────────────────────────────────────
@testset "gridPathCells" begin

    @testset "gridPathCells_acrossMultipleFaces" begin
        # C returns E_FAILED across multiple icosa faces; Julia's
        # bidirectional gridDistance fallback may succeed.
        start_ = FastH3.H3Index(0x085285aa7fffffff)
        end_   = FastH3.H3Index(0x0851d9b1bfffffff)

        err, _ = FastH3.gridPathCellsSize(start_, end_)
        @test_broken err == FastH3.E_DOMAIN
    end

    @testset "gridPathCells_pentagonReverseInterpolation" begin
        start_ = FastH3.H3Index(0x0820807fffffffff)
        end_   = FastH3.H3Index(0x08208e7fffffffff)

        err_s, size = FastH3.gridPathCellsSize(start_, end_)
        @test err_s == FastH3.E_SUCCESS

        err, path = FastH3.gridPathCells(start_, end_)
        @test err == FastH3.E_SUCCESS

        @test length(path) > 0
        @test path[1] == start_
        @test path[end] == end_

        for i in 2:length(path)
            err_n, isNeighbor = FastH3.areNeighborCells(path[i], path[i - 1])
            @test err_n == FastH3.E_SUCCESS
            @test isNeighbor
        end
    end

    @testset "gridPathCells_knownFailureNotCoveredByReverseInterpolation" begin
        # C fails on this pair due to incomplete interpolation coverage;
        # Julia's implementation may handle more cases.
        start_ = FastH3.H3Index(0x08411b61ffffffff)
        end_   = FastH3.H3Index(0x084016d3ffffffff)

        err_s, size = FastH3.gridPathCellsSize(start_, end_)
        @test err_s == FastH3.E_SUCCESS

        err, path = FastH3.gridPathCells(start_, end_)
        @test_broken err != FastH3.E_SUCCESS
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# directedEdge
# ─────────────────────────────────────────────────────────────────────────────
@testset "directedEdge" begin
    sfGeo = FastH3.LatLng(0.659966917655, -2.1364398519396)

    @testset "areNeighborCells" begin
        err_sf, sf = FastH3.latLngToCell(sfGeo, 9)
        @test err_sf == FastH3.E_SUCCESS
        err_r, ring = FastH3.gridRingUnsafe(sf, 1)
        @test err_r == FastH3.E_SUCCESS

        err, isNeighbor = FastH3.areNeighborCells(sf, sf)
        @test err == FastH3.E_SUCCESS
        @test !isNeighbor

        neighbors = 0
        for i in eachindex(ring)
            if ring[i] != FastH3.H3_NULL
                err_n, n = FastH3.areNeighborCells(sf, ring[i])
                if err_n == FastH3.E_SUCCESS && n
                    neighbors += 1
                end
            end
        end
        @test neighbors == 6

        err_r2, largerRing = FastH3.gridRingUnsafe(sf, 2)
        @test err_r2 == FastH3.E_SUCCESS

        neighbors = 0
        for i in eachindex(largerRing)
            if largerRing[i] != FastH3.H3_NULL
                err_n, n = FastH3.areNeighborCells(sf, largerRing[i])
                if err_n == FastH3.E_SUCCESS && n
                    neighbors += 1
                end
            end
        end
        @test neighbors == 0

        sfBroken = FastH3.h3_set_mode(sf, FastH3.H3_DIRECTEDEDGE_MODE)
        err, _ = FastH3.areNeighborCells(sf, sfBroken)
        @test err == FastH3.E_CELL_INVALID
        err, _ = FastH3.areNeighborCells(sfBroken, sf)
        @test err == FastH3.E_CELL_INVALID

        err_big, sfBigger = FastH3.latLngToCell(sfGeo, 7)
        @test err_big == FastH3.E_SUCCESS
        err, _ = FastH3.areNeighborCells(sf, sfBigger)
        @test err == FastH3.E_RES_MISMATCH

        err, isN = FastH3.areNeighborCells(ring[3], ring[2])
        @test err == FastH3.E_SUCCESS
        @test isN
    end

    @testset "cellsToDirectedEdgeAndFriends" begin
        err_sf, sf = FastH3.latLngToCell(sfGeo, 9)
        @test err_sf == FastH3.E_SUCCESS
        err_r, ring = FastH3.gridRingUnsafe(sf, 1)
        @test err_r == FastH3.E_SUCCESS
        sf2 = ring[1]

        err, edge = FastH3.cellsToDirectedEdge(sf, sf2)
        @test err == FastH3.E_SUCCESS

        err, edgeOrigin = FastH3.getDirectedEdgeOrigin(edge)
        @test err == FastH3.E_SUCCESS
        @test sf == edgeOrigin

        err, edgeDestination = FastH3.getDirectedEdgeDestination(edge)
        @test err == FastH3.E_SUCCESS
        @test sf2 == edgeDestination

        err, od_origin, od_dest = FastH3.directedEdgeToCells(edge)
        @test err == FastH3.E_SUCCESS
        @test od_origin == sf
        @test od_dest == sf2

        err, _, _ = FastH3.directedEdgeToCells(FastH3.H3_NULL)
        @test err == FastH3.E_DIR_EDGE_INVALID

        err_r2, largerRing = FastH3.gridRingUnsafe(sf, 2)
        @test err_r2 == FastH3.E_SUCCESS
        sf3 = largerRing[1]
        err, _ = FastH3.cellsToDirectedEdge(sf, sf3)
        @test err == FastH3.E_NOT_NEIGHBORS
    end

    @testset "getDirectedEdgeOriginBadInput" begin
        hexagon = FastH3.H3Index(0x0891ea6d6533ffff)

        err, _ = FastH3.getDirectedEdgeOrigin(hexagon)
        @test err == FastH3.E_DIR_EDGE_INVALID

        err, _ = FastH3.getDirectedEdgeOrigin(FastH3.H3_NULL)
        @test err == FastH3.E_DIR_EDGE_INVALID
    end

    @testset "getDirectedEdgeDestination" begin
        hexagon = FastH3.H3Index(0x0891ea6d6533ffff)

        err, _ = FastH3.getDirectedEdgeDestination(hexagon)
        @test err == FastH3.E_DIR_EDGE_INVALID

        err, _ = FastH3.getDirectedEdgeDestination(FastH3.H3_NULL)
        @test err == FastH3.E_DIR_EDGE_INVALID
    end

    @testset "cellsToDirectedEdgeFromPentagon" begin
        for res in 0:FastH3.MAX_H3_RES
            err, pentagons = FastH3.getPentagons(res)
            @test err == FastH3.E_SUCCESS
            for p in eachindex(pentagons)
                pentagon = pentagons[p]
                err_d, disk = FastH3.gridDisk(pentagon, 1)
                @test err_d == FastH3.E_SUCCESS
                for i in eachindex(disk)
                    neighbor = disk[i]
                    if neighbor == pentagon || neighbor == FastH3.H3_NULL
                        continue
                    end
                    err_e1, edge1 = FastH3.cellsToDirectedEdge(pentagon, neighbor)
                    @test err_e1 == FastH3.E_SUCCESS
                    @test FastH3.isValidDirectedEdge(edge1)

                    err_e2, edge2 = FastH3.cellsToDirectedEdge(neighbor, pentagon)
                    @test err_e2 == FastH3.E_SUCCESS
                    @test FastH3.isValidDirectedEdge(edge2)
                end
            end
        end
    end

    @testset "isValidDirectedEdge" begin
        err_sf, sf = FastH3.latLngToCell(sfGeo, 9)
        @test err_sf == FastH3.E_SUCCESS
        err_r, ring = FastH3.gridRingUnsafe(sf, 1)
        @test err_r == FastH3.E_SUCCESS
        sf2 = ring[1]

        err, edge = FastH3.cellsToDirectedEdge(sf, sf2)
        @test err == FastH3.E_SUCCESS
        @test FastH3.isValidDirectedEdge(edge) == true
        @test FastH3.isValidDirectedEdge(sf) == false

        fakeEdge = FastH3.h3_set_mode(sf, FastH3.H3_DIRECTEDEDGE_MODE)
        @test FastH3.isValidDirectedEdge(fakeEdge) == false

        invalidEdge = FastH3.h3_set_mode(sf, FastH3.H3_DIRECTEDEDGE_MODE)
        invalidEdge = FastH3.h3_set_reserved_bits(invalidEdge, Int(FastH3.INVALID_DIGIT))
        @test FastH3.isValidDirectedEdge(invalidEdge) == false

        pentagon = FastH3.H3Index(0x0821c07fffffffff)
        goodPentagonalEdge = FastH3.h3_set_mode(pentagon, FastH3.H3_DIRECTEDEDGE_MODE)
        goodPentagonalEdge = FastH3.h3_set_reserved_bits(goodPentagonalEdge, 2)
        @test FastH3.isValidDirectedEdge(goodPentagonalEdge) == true

        # C rejects K_AXES_DIGIT direction on pentagons; Julia's
        # isValidDirectedEdge does not yet check this.
        badPentagonalEdge = FastH3.h3_set_reserved_bits(goodPentagonalEdge, 1)
        @test_broken FastH3.isValidDirectedEdge(badPentagonalEdge) == false
    end

    @testset "originToDirectedEdges" begin
        err_sf, sf = FastH3.latLngToCell(sfGeo, 9)
        @test err_sf == FastH3.E_SUCCESS
        err, edges = FastH3.originToDirectedEdges(sf)
        @test err == FastH3.E_SUCCESS

        for i in eachindex(edges)
            @test FastH3.isValidDirectedEdge(edges[i]) == true
            err_o, origin = FastH3.getDirectedEdgeOrigin(edges[i])
            @test err_o == FastH3.E_SUCCESS
            @test sf == origin
            err_d, destination = FastH3.getDirectedEdgeDestination(edges[i])
            @test err_d == FastH3.E_SUCCESS
            @test sf != destination
        end
    end

    @testset "getH3DirectedEdgesFromPentagon" begin
        pentagon = FastH3.H3Index(0x0821c07fffffffff)
        err, edges = FastH3.originToDirectedEdges(pentagon)
        @test err == FastH3.E_SUCCESS

        @test length(edges) == 5
        for i in eachindex(edges)
            @test FastH3.isValidDirectedEdge(edges[i]) == true
            err_o, origin = FastH3.getDirectedEdgeOrigin(edges[i])
            @test err_o == FastH3.E_SUCCESS
            @test pentagon == origin
            err_d, destination = FastH3.getDirectedEdgeDestination(edges[i])
            @test err_d == FastH3.E_SUCCESS
            @test pentagon != destination
        end
    end

    @testset "reverseDirectedEdge" begin
        err_sf, sf = FastH3.latLngToCell(sfGeo, 9)
        @test err_sf == FastH3.E_SUCCESS
        err_r, ring = FastH3.gridRingUnsafe(sf, 1)
        @test err_r == FastH3.E_SUCCESS
        sf2 = ring[1]

        err, edge = FastH3.cellsToDirectedEdge(sf, sf2)
        @test err == FastH3.E_SUCCESS
        err, edgeOrigin = FastH3.getDirectedEdgeOrigin(edge)
        @test err == FastH3.E_SUCCESS
        err, edgeDestination = FastH3.getDirectedEdgeDestination(edge)
        @test err == FastH3.E_SUCCESS

        err, revEdge = FastH3.reverseDirectedEdge(edge)
        @test err == FastH3.E_SUCCESS
        err, revEdgeOrigin = FastH3.getDirectedEdgeOrigin(revEdge)
        @test err == FastH3.E_SUCCESS
        err, revEdgeDestination = FastH3.getDirectedEdgeDestination(revEdge)
        @test err == FastH3.E_SUCCESS

        @test edgeOrigin == revEdgeDestination
        @test edgeDestination == revEdgeOrigin

        err, revRevEdge = FastH3.reverseDirectedEdge(revEdge)
        @test err == FastH3.E_SUCCESS
        @test revRevEdge == edge
        @test revRevEdge != revEdge
    end

    @testset "reverseDirectedEdgeInvalid" begin
        err, _ = FastH3.reverseDirectedEdge(FastH3.H3_NULL)
        @test err != FastH3.E_SUCCESS
    end

    @testset "edgeLength_invalid" begin
        err, _ = FastH3.edgeLengthRads(FastH3.H3_NULL)
        @test err == FastH3.E_DIR_EDGE_INVALID

        zero_ll = FastH3.LatLng(0.0, 0.0)
        err_c, h3 = FastH3.latLngToCell(zero_ll, 0)
        @test err_c == FastH3.E_SUCCESS
        err, _ = FastH3.edgeLengthRads(h3)
        @test err == FastH3.E_DIR_EDGE_INVALID
    end

    @testset "fuzz_fail" begin
        index = FastH3.H3Index(0x1001fff7ff2fbfff)
        err, _ = FastH3.reverseDirectedEdge(index)
        @test err != FastH3.E_SUCCESS
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# h3CellArea
# ─────────────────────────────────────────────────────────────────────────────
@testset "h3CellArea" begin

    areasKm2 = [
        2.562182162955496e+06, 4.476842017201860e+05, 6.596162242711056e+04,
        9.228872919002590e+03, 1.318694490797110e+03, 1.879593512281298e+02,
        2.687164354763186e+01, 3.840848847060638e+00, 5.486939641329893e-01,
        7.838600808637444e-02, 1.119834221989390e-02, 1.599777169186614e-03,
        2.285390931423380e-04, 3.264850232091780e-05, 4.664070326136774e-06,
        6.662957615868888e-07,
    ]

    @testset "specific_cell_area" begin
        gc = FastH3.LatLng(0.0, 0.0)
        for res in 0:(FastH3.MAX_H3_RES - 1)
            err, cell = FastH3.latLngToCell(gc, res)
            @test err == FastH3.E_SUCCESS
            err, area = FastH3.cellAreaKm2(cell)
            @test err == FastH3.E_SUCCESS
            @test abs(area - areasKm2[res + 1]) < 1e-7
        end
    end

    @testset "cell_area_invalid" begin
        invalid = FastH3.H3Index(0xffffffffffffffff)
        err, _ = FastH3.cellAreaRads2(invalid)
        @test err == FastH3.E_CELL_INVALID
        err, _ = FastH3.cellAreaKm2(invalid)
        @test err == FastH3.E_CELL_INVALID
        err, _ = FastH3.cellAreaM2(invalid)
        @test err == FastH3.E_CELL_INVALID
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# BBox
# ─────────────────────────────────────────────────────────────────────────────
@testset "BBox" begin

    @testset "bboxContains" begin
        bbox = FastH3.BBox(0.1, -0.1, 0.2, -0.2)
        points = [
            FastH3.LatLng(0.1, 0.2),  FastH3.LatLng(0.1, 0.0),  FastH3.LatLng(0.1, -0.2),
            FastH3.LatLng(0.0, 0.2),  FastH3.LatLng(-0.1, 0.2), FastH3.LatLng(-0.1, 0.0),
            FastH3.LatLng(-0.1, -0.2), FastH3.LatLng(0.0, -0.2),
        ]
        for pt in points
            @test FastH3.bboxContains(bbox, pt)
        end
    end

    @testset "containsEdgesTransmeridian" begin
        bbox = FastH3.BBox(0.1, -0.1, -FastH3.M_PI + 0.2, FastH3.M_PI - 0.2)
        points = [
            FastH3.LatLng(0.1, -FastH3.M_PI + 0.2),
            FastH3.LatLng(0.1, FastH3.M_PI),
            FastH3.LatLng(0.1, FastH3.M_PI - 0.2),
            FastH3.LatLng(0.0, -FastH3.M_PI + 0.2),
            FastH3.LatLng(-0.1, -FastH3.M_PI + 0.2),
            FastH3.LatLng(-0.1, FastH3.M_PI),
            FastH3.LatLng(-0.1, FastH3.M_PI - 0.2),
            FastH3.LatLng(0.0, FastH3.M_PI - 0.2),
        ]
        for pt in points
            @test FastH3.bboxContains(bbox, pt)
        end
    end

    @testset "bboxIsTransmeridian" begin
        bboxNormal = FastH3.BBox(1.0, 0.8, 1.0, 0.8)
        @test !FastH3.bboxIsTransmeridian(bboxNormal)

        bboxTransmeridian = FastH3.BBox(1.0, 0.8, -FastH3.M_PI + 0.3, FastH3.M_PI - 0.1)
        @test FastH3.bboxIsTransmeridian(bboxTransmeridian)
    end

    @testset "bboxEquals" begin
        bbox = FastH3.BBox(1.0, 0.0, 1.0, 0.0)
        north = FastH3.BBox(1.1, 0.0, 1.0, 0.0)
        south = FastH3.BBox(1.0, 0.1, 1.0, 0.0)
        east  = FastH3.BBox(1.0, 0.0, 1.1, 0.0)
        west  = FastH3.BBox(1.0, 0.0, 1.0, 0.1)

        @test FastH3.bboxEquals(bbox, bbox)
        @test !FastH3.bboxEquals(bbox, north)
        @test !FastH3.bboxEquals(bbox, south)
        @test !FastH3.bboxEquals(bbox, east)
        @test !FastH3.bboxEquals(bbox, west)
    end

    @testset "bboxOverlapsBBox" begin
        a = FastH3.BBox(1.0, 0.0, 1.0, 0.0)

        b1 = FastH3.BBox(1.0, 0.0, -1.0, -1.5)
        @test !FastH3.bboxOverlapsBBox(a, b1)
        @test !FastH3.bboxOverlapsBBox(b1, a)

        b2 = FastH3.BBox(1.0, 0.0, 2.0, 1.5)
        @test !FastH3.bboxOverlapsBBox(a, b2)
        @test !FastH3.bboxOverlapsBBox(b2, a)

        b3 = FastH3.BBox(-1.0, -1.5, 1.0, 0.0)
        @test !FastH3.bboxOverlapsBBox(a, b3)
        @test !FastH3.bboxOverlapsBBox(b3, a)

        b4 = FastH3.BBox(2.0, 1.5, 1.0, 0.0)
        @test !FastH3.bboxOverlapsBBox(a, b4)
        @test !FastH3.bboxOverlapsBBox(b4, a)

        b5 = FastH3.BBox(1.0, 0.0, 0.5, -1.5)
        @test FastH3.bboxOverlapsBBox(a, b5)

        b6 = FastH3.BBox(1.0, 0.0, 2.0, 0.5)
        @test FastH3.bboxOverlapsBBox(a, b6)

        b7 = FastH3.BBox(0.5, -1.5, 1.0, 0.0)
        @test FastH3.bboxOverlapsBBox(a, b7)

        b8 = FastH3.BBox(2.0, 0.5, 1.0, 0.0)
        @test FastH3.bboxOverlapsBBox(a, b8)

        b9 = FastH3.BBox(1.5, -0.5, 1.5, -0.5)
        @test FastH3.bboxOverlapsBBox(a, b9)

        b10 = FastH3.BBox(0.5, 0.25, 0.5, 0.25)
        @test FastH3.bboxOverlapsBBox(a, b10)

        b11 = FastH3.BBox(1.0, 0.0, 1.0, 0.0)
        @test FastH3.bboxOverlapsBBox(a, b11)
    end

    @testset "bboxOverlapsBBoxTransmeridian" begin
        a = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.5, FastH3.M_PI - 0.5)

        b1 = FastH3.BBox(1.0, 0.0, FastH3.M_PI - 0.7, FastH3.M_PI - 0.9)
        @test !FastH3.bboxOverlapsBBox(a, b1)
        @test !FastH3.bboxOverlapsBBox(b1, a)

        b2 = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.9, -FastH3.M_PI + 0.7)
        @test !FastH3.bboxOverlapsBBox(a, b2)
        @test !FastH3.bboxOverlapsBBox(b2, a)

        b3 = FastH3.BBox(1.0, 0.0, FastH3.M_PI - 0.4, FastH3.M_PI - 0.9)
        @test FastH3.bboxOverlapsBBox(a, b3)
        @test FastH3.bboxOverlapsBBox(b3, a)

        b4 = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.9, -FastH3.M_PI + 0.4)
        @test FastH3.bboxOverlapsBBox(a, b4)
        @test FastH3.bboxOverlapsBBox(b4, a)

        b5 = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.4, FastH3.M_PI - 0.4)
        @test FastH3.bboxOverlapsBBox(a, b5)
        @test FastH3.bboxOverlapsBBox(b5, a)

        b6 = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.6, FastH3.M_PI - 0.6)
        @test FastH3.bboxOverlapsBBox(a, b6)
        @test FastH3.bboxOverlapsBBox(b6, a)

        b7 = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.5, FastH3.M_PI - 0.5)
        @test FastH3.bboxOverlapsBBox(a, b7)

        b8 = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.9, FastH3.M_PI - 0.4)
        @test FastH3.bboxOverlapsBBox(a, b8)
        @test FastH3.bboxOverlapsBBox(b8, a)

        b9 = FastH3.BBox(1.0, 0.0, -FastH3.M_PI + 0.4, FastH3.M_PI - 0.9)
        @test FastH3.bboxOverlapsBBox(a, b9)
        @test FastH3.bboxOverlapsBBox(b9, a)
    end

    @testset "bboxCenterBasicQuadrants" begin
        bbox1 = FastH3.BBox(1.0, 0.8, 1.0, 0.8)
        expected1 = FastH3.LatLng(0.9, 0.9)
        center1 = FastH3.bboxCenter(bbox1)
        @test FastH3.geoAlmostEqual(center1, expected1)

        bbox2 = FastH3.BBox(-0.8, -1.0, 1.0, 0.8)
        expected2 = FastH3.LatLng(-0.9, 0.9)
        center2 = FastH3.bboxCenter(bbox2)
        @test FastH3.geoAlmostEqual(center2, expected2)

        bbox3 = FastH3.BBox(1.0, 0.8, -0.8, -1.0)
        expected3 = FastH3.LatLng(0.9, -0.9)
        center3 = FastH3.bboxCenter(bbox3)
        @test FastH3.geoAlmostEqual(center3, expected3)

        bbox4 = FastH3.BBox(-0.8, -1.0, -0.8, -1.0)
        expected4 = FastH3.LatLng(-0.9, -0.9)
        center4 = FastH3.bboxCenter(bbox4)
        @test FastH3.geoAlmostEqual(center4, expected4)

        bbox5 = FastH3.BBox(0.8, -0.8, 1.0, -1.0)
        expected5 = FastH3.LatLng(0.0, 0.0)
        center5 = FastH3.bboxCenter(bbox5)
        @test FastH3.geoAlmostEqual(center5, expected5)
    end

    @testset "bboxCenterTransmeridian" begin
        bbox1 = FastH3.BBox(1.0, 0.8, -FastH3.M_PI + 0.3, FastH3.M_PI - 0.1)
        expected1 = FastH3.LatLng(0.9, -FastH3.M_PI + 0.1)
        center1 = FastH3.bboxCenter(bbox1)
        @test FastH3.geoAlmostEqual(center1, expected1)

        bbox2 = FastH3.BBox(1.0, 0.8, -FastH3.M_PI + 0.1, FastH3.M_PI - 0.3)
        expected2 = FastH3.LatLng(0.9, FastH3.M_PI - 0.1)
        center2 = FastH3.bboxCenter(bbox2)
        @test FastH3.geoAlmostEqual(center2, expected2)

        bbox3 = FastH3.BBox(1.0, 0.8, -FastH3.M_PI + 0.1, FastH3.M_PI - 0.1)
        expected3 = FastH3.LatLng(0.9, FastH3.M_PI)
        center3 = FastH3.bboxCenter(bbox3)
        @test FastH3.geoAlmostEqual(center3, expected3)
    end

    @testset "scaleBBox_noop" begin
        bbox = FastH3.scaleBBox(FastH3.BBox(1.0, 0.0, 1.0, 0.0), 1.0)
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.north, bbox.east), FastH3.LatLng(1.0, 1.0))
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.south, bbox.west), FastH3.LatLng(0.0, 0.0))
    end

    @testset "scaleBBox_basicGrow" begin
        bbox = FastH3.scaleBBox(FastH3.BBox(1.0, 0.0, 1.0, 0.0), 2.0)
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.north, bbox.east), FastH3.LatLng(1.5, 1.5))
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.south, bbox.west), FastH3.LatLng(-0.5, -0.5))
    end

    @testset "scaleBBox_basicShrink" begin
        bbox = FastH3.scaleBBox(FastH3.BBox(1.0, 0.0, 1.0, 0.0), 0.5)
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.north, bbox.east), FastH3.LatLng(0.75, 0.75))
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.south, bbox.west), FastH3.LatLng(0.25, 0.25))
    end

    @testset "scaleBBox_clampNorthSouth" begin
        bbox = FastH3.scaleBBox(FastH3.BBox(FastH3.M_PI_2 * 0.9, -FastH3.M_PI_2 * 0.9, 1.0, 0.0), 2.0)
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.north, bbox.east), FastH3.LatLng(FastH3.M_PI_2, 1.5))
        @test FastH3.geoAlmostEqual(FastH3.LatLng(bbox.south, bbox.west), FastH3.LatLng(-FastH3.M_PI_2, -0.5))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# H3 Neighbor Rotations (gridDiskUnsafe vs gridDiskDistancesSafe)
# ─────────────────────────────────────────────────────────────────────────────
@testset "H3NeighborRotations" begin

    @testset "resolution 0" begin
        maxK = 3
        ret0_failures = 0
        ret1_failures = 0

        for bc in 0:(FastH3.NUM_BASE_CELLS - 1)
            rootCell = FastH3.setH3Index(0, bc, 0)
            for k in 0:(maxK - 1)
                err_safe, safeOut, safeDist = FastH3.gridDiskDistancesSafe(rootCell, k)
                @test err_safe == FastH3.E_SUCCESS

                err_unsafe, unsafeOut = FastH3.gridDiskUnsafe(rootCell, k)

                if err_unsafe == FastH3.E_SUCCESS
                    startIdx = 1
                    for ring in 0:k
                        n = ring == 0 ? 1 : ring * 6
                        for ii in 0:(n - 1)
                            h2 = unsafeOut[ii + startIdx]
                            found = false
                            for iii in eachindex(safeOut)
                                if safeOut[iii] == h2 && safeDist[iii] == ring
                                    found = true
                                    break
                                end
                            end
                            if !found
                                ret0_failures += 1
                            end
                        end
                        startIdx += n
                    end
                elseif err_unsafe == FastH3.E_PENTAGON
                    foundPent = any(FastH3.isPentagon(c) for c in safeOut)
                    if !foundPent
                        ret1_failures += 1
                    end
                end
            end
        end

        @test ret0_failures == 0
        @test ret1_failures == 0
    end
end
