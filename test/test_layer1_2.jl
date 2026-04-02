# Tests are included from runtests.jl which loads H3X

# ============================================================================
# testCoordIjkInternal
# ============================================================================
@testset "coordIjkInternal" begin
    @testset "_unitIjkToDigit" begin
        zero = H3X.CoordIJK(0, 0, 0)
        i = H3X.CoordIJK(1, 0, 0)
        outOfRange = H3X.CoordIJK(2, 0, 0)
        unnormalizedZero = H3X.CoordIJK(2, 2, 2)

        @test H3X._unitIjkToDigit(zero) == H3X.CENTER_DIGIT
        @test H3X._unitIjkToDigit(i) == H3X.I_AXES_DIGIT
        @test H3X._unitIjkToDigit(outOfRange) == H3X.INVALID_DIGIT
        @test H3X._unitIjkToDigit(unnormalizedZero) == H3X.CENTER_DIGIT
    end

    @testset "_neighbor" begin
        ijk = H3X.CoordIJK(0, 0, 0)

        zero = H3X.CoordIJK(0, 0, 0)
        i_coord = H3X.CoordIJK(1, 0, 0)

        ijk = H3X._neighbor(ijk, H3X.CENTER_DIGIT)
        @test H3X._ijkMatches(ijk, zero)
        ijk = H3X._neighbor(ijk, H3X.I_AXES_DIGIT)
        @test H3X._ijkMatches(ijk, i_coord)
        ijk = H3X._neighbor(ijk, H3X.INVALID_DIGIT)
        @test H3X._ijkMatches(ijk, i_coord)
    end

    @testset "_upAp7Checked" begin
        @test first(H3X._upAp7Checked(H3X.CoordIJK(0, 0, 0))) == H3X.E_SUCCESS
        @test first(H3X._upAp7Checked(H3X.CoordIJK(typemax(Int32), Int32(0), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7Checked(H3X.CoordIJK(typemax(Int32) ÷ Int32(2), Int32(0), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7Checked(H3X.CoordIJK(Int32(0), typemax(Int32), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7Checked(H3X.CoordIJK(typemax(Int32) ÷ Int32(3), Int32(-2), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7Checked(H3X.CoordIJK(typemax(Int32) ÷ Int32(3), typemax(Int32) ÷ Int32(2), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7Checked(H3X.CoordIJK(Int32(-1), Int32(0), Int32(0)))) == H3X.E_SUCCESS
    end

    @testset "_upAp7rChecked" begin
        @test first(H3X._upAp7rChecked(H3X.CoordIJK(0, 0, 0))) == H3X.E_SUCCESS
        @test first(H3X._upAp7rChecked(H3X.CoordIJK(typemax(Int32), Int32(0), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7rChecked(H3X.CoordIJK(Int32(0), typemax(Int32), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7rChecked(H3X.CoordIJK(Int32(0), typemax(Int32) ÷ Int32(2), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7rChecked(H3X.CoordIJK(typemax(Int32) ÷ Int32(2), typemax(Int32) ÷ Int32(3), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7rChecked(H3X.CoordIJK(Int32(-2), typemax(Int32) ÷ Int32(3), Int32(0)))) == H3X.E_FAILED
        @test first(H3X._upAp7rChecked(H3X.CoordIJK(Int32(-1), Int32(0), Int32(0)))) == H3X.E_SUCCESS
    end
end

# ============================================================================
# testBaseCells
# ============================================================================
@testset "baseCells" begin
    @testset "getRes0Cells" begin
        (err, indexes) = H3X.getRes0Cells()
        @test err == H3X.E_SUCCESS
        @test indexes[1] == H3X.H3Index(0x08001fffffffffff)
        @test indexes[122] == H3X.H3Index(0x080f3fffffffffff)
    end
end

# ============================================================================
# testBaseCellsInternal
# ============================================================================
@testset "baseCellsInternal" begin
    @testset "baseCellToCCWrot60" begin
        @test H3X._baseCellToCCWrot60(16, 0) == 0
        @test H3X._baseCellToCCWrot60(32, 0) == 3
        @test H3X._baseCellToCCWrot60(7, 3) == 1
    end

    @testset "baseCellToCCWrot60_invalid" begin
        @test H3X._baseCellToCCWrot60(16, 42) == H3X.INVALID_ROTATIONS
        @test H3X._baseCellToCCWrot60(16, -1) == H3X.INVALID_ROTATIONS
        @test H3X._baseCellToCCWrot60(1, 0) == H3X.INVALID_ROTATIONS
    end

    @testset "isBaseCellPentagon_invalid" begin
        @test H3X._isBaseCellPentagon(-1) == false
    end
end

# ============================================================================
# testH3Index
# ============================================================================
@testset "h3Index" begin
    @testset "latLngToCellExtremeCoordinates" begin
        g = H3X.LatLng(0.0, 1e45)
        (err, h) = H3X.latLngToCell(g, 14)
        @test err == H3X.E_SUCCESS

        g2 = H3X.LatLng(1e46, 1e45)
        (err2, h2) = H3X.latLngToCell(g2, 15)
        @test err2 == H3X.E_SUCCESS

        g4 = H3X.setGeoDegs(2.0, -3e39)
        (err4, h4) = H3X.latLngToCell(g4, 0)
        @test err4 == H3X.E_SUCCESS
    end

    @testset "isValidCellAtResolution" begin
        for i in 0:H3X.MAX_H3_RES
            g = H3X.LatLng(0.0, 0.0)
            (err, h3) = H3X.latLngToCell(g, i)
            @test err == H3X.E_SUCCESS
            @test H3X.isValidCell(h3)
        end
    end

    @testset "isValidCellDigits" begin
        g = H3X.LatLng(0.0, 0.0)
        (err, h3) = H3X.latLngToCell(g, 1)
        @test err == H3X.E_SUCCESS
        h3 = xor(h3, UInt64(1))
        @test !H3X.isValidCell(h3)
    end

    @testset "isValidCellBaseCell" begin
        for i in 0:(H3X.NUM_BASE_CELLS - 1)
            h = H3X.H3_INIT
            h = H3X.h3_set_mode(h, H3X.H3_CELL_MODE)
            h = H3X.h3_set_base_cell(h, i)
            @test H3X.isValidCell(h)
            @test H3X.getBaseCellNumber(h) == i
        end
    end

    @testset "isValidCellBaseCellInvalid" begin
        hWrongBaseCell = H3X.H3_INIT
        hWrongBaseCell = H3X.h3_set_mode(hWrongBaseCell, H3X.H3_CELL_MODE)
        hWrongBaseCell = H3X.h3_set_base_cell(hWrongBaseCell, H3X.NUM_BASE_CELLS)
        @test !H3X.isValidCell(hWrongBaseCell)
    end

    @testset "isValidCellWithMode" begin
        for i in 0:0xf
            h = H3X.H3_INIT
            h = H3X.h3_set_mode(h, i)
            if i == H3X.H3_CELL_MODE
                @test H3X.isValidCell(h)
            else
                @test !H3X.isValidCell(h)
            end
        end
    end

    @testset "isValidCellReservedBits" begin
        for i in 0:7
            h = H3X.H3_INIT
            h = H3X.h3_set_mode(h, H3X.H3_CELL_MODE)
            h = H3X.h3_set_reserved_bits(h, i)
            if i == 0
                @test H3X.isValidCell(h)
            else
                @test !H3X.isValidCell(h)
            end
        end
    end

    @testset "isValidCellHighBit" begin
        h = H3X.H3_INIT
        h = H3X.h3_set_mode(h, H3X.H3_CELL_MODE)
        h = h | H3X.H3_HIGH_BIT_MASK
        @test !H3X.isValidCell(h)
    end

    @testset "h3BadDigitInvalid" begin
        h = H3X.H3_INIT
        h = H3X.h3_set_mode(h, H3X.H3_CELL_MODE)
        h = H3X.h3_set_resolution(h, 1)
        @test !H3X.isValidCell(h)
    end

    @testset "h3DeletedSubsequenceInvalid" begin
        h = H3X.setH3Index(1, 4, H3X.K_AXES_DIGIT)
        @test !H3X.isValidCell(h)
    end

    @testset "moreDeletedSubsequenceInvalid" begin
        p = H3X.H3Index(0x080c3fffffffffff)  # res 0 pentagon
        for res in 1:15
            (err, h) = H3X.cellToCenterChild(p, res)
            @test err == H3X.E_SUCCESS
            @test H3X.isValidCell(h)
            for d in 0:6
                h2 = H3X.h3_set_index_digit(h, res, d)
                if d == 1
                    @test !H3X.isValidCell(h2)
                else
                    @test H3X.isValidCell(h2)
                end
            end
        end
    end

    @testset "h3ToString" begin
        @test H3X.h3ToString(UInt64(0xcafe)) == "cafe"
        @test H3X.h3ToString(typemax(UInt64)) == "ffffffffffffffff"
    end

    @testset "stringToH3" begin
        (err, h3) = H3X.stringToH3("")
        @test err == H3X.E_FAILED

        (err2, h3_2) = H3X.stringToH3("**")
        @test err2 == H3X.E_FAILED

        (err3, h3_3) = H3X.stringToH3("ffffffffffffffff")
        @test err3 == H3X.E_SUCCESS
        @test h3_3 == 0xffffffffffffffff
    end

    @testset "setH3Index" begin
        h = H3X.setH3Index(5, 12, 1)
        @test H3X.h3_get_resolution(h) == 5
        @test H3X.h3_get_base_cell(h) == 12
        @test H3X.h3_get_mode(h) == H3X.H3_CELL_MODE
        for i in 1:5
            @test H3X.h3_get_index_digit(h, i) == H3X.Direction(1)
        end
        for i in 6:H3X.MAX_H3_RES
            @test H3X.h3_get_index_digit(h, i) == H3X.Direction(7)
        end
        @test h == H3X.H3Index(0x085184927fffffff)
    end

    @testset "isResClassIII" begin
        coord = H3X.LatLng(0.0, 0.0)
        for i in 0:H3X.MAX_H3_RES
            (err, h) = H3X.latLngToCell(coord, i)
            @test err == H3X.E_SUCCESS
            @test H3X.isResClassIII(h) == H3X.isResolutionClassIII(i)
        end
    end
end

# ============================================================================
# testH3IndexInternal
# ============================================================================
@testset "h3IndexInternal" begin
    @testset "faceIjkToH3ExtremeCoordinates" begin
        fijk0I = H3X.FaceIJK(0, H3X.CoordIJK(3, 0, 0))
        @test H3X._faceIjkToH3(fijk0I, 0) == H3X.H3Index(0)
        fijk0J = H3X.FaceIJK(1, H3X.CoordIJK(0, 4, 0))
        @test H3X._faceIjkToH3(fijk0J, 0) == H3X.H3Index(0)
        fijk0K = H3X.FaceIJK(2, H3X.CoordIJK(2, 0, 5))
        @test H3X._faceIjkToH3(fijk0K, 0) == H3X.H3Index(0)

        fijk1I = H3X.FaceIJK(3, H3X.CoordIJK(6, 0, 0))
        @test H3X._faceIjkToH3(fijk1I, 1) == H3X.H3Index(0)
        fijk1J = H3X.FaceIJK(4, H3X.CoordIJK(0, 7, 1))
        @test H3X._faceIjkToH3(fijk1J, 1) == H3X.H3Index(0)
        fijk1K = H3X.FaceIJK(5, H3X.CoordIJK(2, 0, 8))
        @test H3X._faceIjkToH3(fijk1K, 1) == H3X.H3Index(0)

        fijk2I = H3X.FaceIJK(6, H3X.CoordIJK(18, 0, 0))
        @test H3X._faceIjkToH3(fijk2I, 2) == H3X.H3Index(0)
        fijk2J = H3X.FaceIJK(7, H3X.CoordIJK(0, 19, 1))
        @test H3X._faceIjkToH3(fijk2J, 2) == H3X.H3Index(0)
        fijk2K = H3X.FaceIJK(8, H3X.CoordIJK(2, 0, 20))
        @test H3X._faceIjkToH3(fijk2K, 2) == H3X.H3Index(0)
    end
end

# ============================================================================
# testH3Api
# ============================================================================
@testset "h3Api" begin
    @testset "latLngToCell_res" begin
        anywhere = H3X.LatLng(0.0, 0.0)
        (err1, _) = H3X.latLngToCell(anywhere, -1)
        @test err1 == H3X.E_RES_DOMAIN
        (err2, _) = H3X.latLngToCell(anywhere, 16)
        @test err2 == H3X.E_RES_DOMAIN
    end

    @testset "latLngToCell_coord" begin
        invalidLat = H3X.LatLng(NaN, 0.0)
        invalidLng = H3X.LatLng(0.0, NaN)
        invalidLatLng = H3X.LatLng(Inf, -Inf)

        (err1, _) = H3X.latLngToCell(invalidLat, 1)
        @test err1 == H3X.E_LATLNG_DOMAIN
        (err2, _) = H3X.latLngToCell(invalidLng, 1)
        @test err2 == H3X.E_LATLNG_DOMAIN
        (err3, _) = H3X.latLngToCell(invalidLatLng, 1)
        @test err3 == H3X.E_LATLNG_DOMAIN
    end

    @testset "cellToBoundary_classIIIEdgeVertex" begin
        hexes = H3X.H3Index[
            0x0894cc5349b7ffff, 0x0894cc534d97ffff, 0x0894cc53682bffff,
            0x0894cc536b17ffff, 0x0894cc53688bffff, 0x0894cead92cbffff,
            0x0894cc536537ffff, 0x0894cc5acbabffff, 0x0894cc536597ffff,
        ]
        for hex in hexes
            (err, b) = H3X.cellToBoundary(hex)
            @test err == H3X.E_SUCCESS
            @test b.numVerts == 7
        end
    end

    @testset "cellToBoundary_classIIIEdgeVertex_exact" begin
        (err, h3) = H3X.stringToH3("894cc536537ffff")
        @test err == H3X.E_SUCCESS

        (err2, boundary) = H3X.cellToBoundary(h3)
        @test err2 == H3X.E_SUCCESS
        @test boundary.numVerts == 7

        expected_verts = [
            H3X.setGeoDegs(18.043333154, -66.27836523500002),
            H3X.setGeoDegs(18.042238363, -66.27929062800001),
            H3X.setGeoDegs(18.040818259, -66.27854193899998),
            H3X.setGeoDegs(18.040492975, -66.27686786700002),
            H3X.setGeoDegs(18.041040385, -66.27640518300001),
            H3X.setGeoDegs(18.041757122, -66.27596711500001),
            H3X.setGeoDegs(18.043007860, -66.27669118199998),
        ]
        epsilon = 0.000001 * H3X.M_PI_180
        for v in 1:7
            @test abs(boundary.verts[v].lat - expected_verts[v].lat) < epsilon
            @test abs(boundary.verts[v].lng - expected_verts[v].lng) < epsilon
        end
    end

    @testset "cellToBoundary_coslngConstrain" begin
        h3 = H3X.H3Index(0x087dc6d364ffffff)
        (err, boundary) = H3X.cellToBoundary(h3)
        @test err == H3X.E_SUCCESS
        @test boundary.numVerts == 6

        expected_verts = [
            H3X.setGeoDegs(-52.0130533678236091, -34.6232931343713091),
            H3X.setGeoDegs(-52.0041156384652012, -34.6096733160584549),
            H3X.setGeoDegs(-51.9929610229502472, -34.6165157145896387),
            H3X.setGeoDegs(-51.9907410568096608, -34.6369680004259877),
            H3X.setGeoDegs(-51.9996738734672377, -34.6505896528323660),
            H3X.setGeoDegs(-52.0108315681413629, -34.6437571897165668),
        ]
        epsilon = 0.000001 * H3X.M_PI_180
        for v in 1:6
            @test abs(boundary.verts[v].lat - expected_verts[v].lat) < epsilon
            @test abs(boundary.verts[v].lng - expected_verts[v].lng) < epsilon
        end
    end

    @testset "cellToBoundary_failed" begin
        h = H3X.H3Index(0x087dc6d364ffffff)
        h = H3X.h3_set_base_cell(h, H3X.NUM_BASE_CELLS + 1)
        (err, _) = H3X.cellToBoundary(h)
        @test err == H3X.E_CELL_INVALID
    end

    @testset "cellToLatLngInvalid" begin
        (err, _) = H3X.cellToLatLng(H3X.H3Index(0x7fffffffffffffff))
        @test err == H3X.E_CELL_INVALID
    end
end

# ============================================================================
# testConstructCell
# ============================================================================
@testset "constructCell" begin
    @testset "tableOfTests" begin
        test_cases = [
            # (expected_result, res, bc, digits)
            # valid cell constructions
            (H3X.H3Index(0x08001fffffffffff), 0, 0, Int[]),
            (H3X.H3Index(0x08003fffffffffff), 0, 1, Int[]),
            (H3X.H3Index(0x080f3fffffffffff), 0, 121, Int[]),
            (H3X.H3Index(0x0839253fffffffff), 3, 73, Int[1, 2, 3]),
            (H3X.H3Index(0x0821f67fffffffff), 2, 15, Int[5, 4]),
            (H3X.H3Index(0x08155bffffffffff), 1, 42, Int[6]),
            (H3X.H3Index(0x08f754e64992d6d8), 15, 58, Int[5, 1, 6, 3, 1, 1, 1, 4, 4, 5, 5, 3, 3, 3, 0]),

            # resolution tests
            (UInt64(H3X.E_RES_DOMAIN), 16, 0, Int[]),
            (UInt64(H3X.E_RES_DOMAIN), 18, 0, Int[]),
            (UInt64(H3X.E_RES_DOMAIN), -1, 0, Int[]),
            (H3X.H3Index(0x08001fffffffffff), 0, 0, Int[]),

            # base cell tests
            (UInt64(H3X.E_BASE_CELL_DOMAIN), 0, 122, Int[]),
            (UInt64(H3X.E_BASE_CELL_DOMAIN), 0, -1, Int[]),
            (UInt64(H3X.E_BASE_CELL_DOMAIN), 0, 259, Int[]),
            (UInt64(H3X.E_BASE_CELL_DOMAIN), 2, 122, Int[1, 0]),

            # digit tests
            (UInt64(H3X.E_DIGIT_DOMAIN), 1, 40, Int[-1]),
            (UInt64(H3X.E_DIGIT_DOMAIN), 1, 40, Int[7]),
            (UInt64(H3X.E_DIGIT_DOMAIN), 1, 40, Int[8]),
            (UInt64(H3X.E_DIGIT_DOMAIN), 1, 40, Int[17]),

            # deleted subsequence tests (bc=4 is pentagon)
            (H3X.H3Index(0x0830800fffffffff), 3, 4, Int[0, 0, 0]),
            (UInt64(H3X.E_DELETED_DIGIT), 3, 4, Int[0, 0, 1]),
            (H3X.H3Index(0x0830802fffffffff), 3, 4, Int[0, 0, 2]),

            # bc=5 is NOT pentagon
            (H3X.H3Index(0x0830a00fffffffff), 3, 5, Int[0, 0, 0]),
            (H3X.H3Index(0x0830a01fffffffff), 3, 5, Int[0, 0, 1]),
            (H3X.H3Index(0x0830a02fffffffff), 3, 5, Int[0, 0, 2]),
        ]

        for (x, res, bc, digits) in test_cases
            valid_tcx = H3X.isValidCell(x)
            (err, h) = H3X.constructCell(res, bc, digits)

            got_expected_valid_cell = valid_tcx && (err == H3X.E_SUCCESS) && (x == h)
            got_expected_error = !valid_tcx && (x == UInt64(err))

            @test got_expected_valid_cell || got_expected_error
        end
    end

    @testset "roundtrip" begin
        for res in 0:4
            all_passed = true
            (err, cells0) = H3X.getRes0Cells()
            @test err == H3X.E_SUCCESS
            for base_h in cells0
                if res == 0
                    h = base_h
                    r = H3X.getResolution(h)
                    bc = H3X.getBaseCellNumber(h)
                    digits = Int[]
                    for rr in 1:r
                        (e, d) = H3X.getIndexDigit(h, rr)
                        e != H3X.E_SUCCESS && (all_passed = false; continue)
                        push!(digits, d)
                    end
                    (e2, out) = H3X.constructCell(r, bc, digits)
                    if e2 != H3X.E_SUCCESS || out != h
                        all_passed = false
                    end
                else
                    (err2, children) = H3X.cellToChildren(base_h, res)
                    if err2 != H3X.E_SUCCESS
                        all_passed = false
                        continue
                    end
                    for h in children
                        r = H3X.getResolution(h)
                        bc = H3X.getBaseCellNumber(h)
                        digits = Int[]
                        for rr in 1:r
                            (e, d) = H3X.getIndexDigit(h, rr)
                            e != H3X.E_SUCCESS && (all_passed = false; continue)
                            push!(digits, d)
                        end
                        (e2, out) = H3X.constructCell(r, bc, digits)
                        if e2 != H3X.E_SUCCESS || out != h
                            all_passed = false
                        end
                    end
                end
            end
            @test all_passed
        end
    end
end

# ============================================================================
# testCellToParent
# ============================================================================
@testset "cellToParent" begin
    sf = H3X.LatLng(0.659966917655, 2 * 3.14159 - 2.1364398519396)

    @testset "ancestorsForEachRes" begin
        for res in 1:14
            for step in 0:(res - 1)
                (err, child) = H3X.latLngToCell(sf, res)
                @test err == H3X.E_SUCCESS
                (err2, parent) = H3X.cellToParent(child, res - step)
                @test err2 == H3X.E_SUCCESS
                (err3, comparisonParent) = H3X.latLngToCell(sf, res - step)
                @test err3 == H3X.E_SUCCESS
                @test parent == comparisonParent
            end
        end
    end

    @testset "invalidInputs" begin
        (err, child) = H3X.latLngToCell(sf, 5)
        @test err == H3X.E_SUCCESS

        (err2, _) = H3X.cellToParent(child, 6)
        @test err2 == H3X.E_RES_MISMATCH
        (err3, _) = H3X.cellToParent(child, -1)
        @test err3 == H3X.E_RES_DOMAIN
        (err4, _) = H3X.cellToParent(child, 15)
        @test err4 == H3X.E_RES_MISMATCH
        (err5, _) = H3X.cellToParent(child, 16)
        @test err5 == H3X.E_RES_DOMAIN
    end
end

# ============================================================================
# testCellToChildren
# ============================================================================
@testset "cellToChildren" begin
    @testset "oneResStep" begin
        h = H3X.H3Index(0x088283080ddfffff)
        (err, children) = H3X.cellToChildren(h, 9)
        @test err == H3X.E_SUCCESS
        expected = H3X.H3Index[
            0x089283080dc3ffff, 0x089283080dc7ffff,
            0x089283080dcbffff, 0x089283080dcfffff,
            0x089283080dd3ffff, 0x089283080dd7ffff,
            0x089283080ddbffff,
        ]
        @test Set(children) == Set(expected)
    end

    @testset "multipleResSteps" begin
        h = H3X.H3Index(0x088283080ddfffff)
        (err, children) = H3X.cellToChildren(h, 10)
        @test err == H3X.E_SUCCESS
        expected = H3X.H3Index[
            0x08a283080dd27fff, 0x08a283080dd37fff, 0x08a283080dc47fff,
            0x08a283080dcdffff, 0x08a283080dc5ffff, 0x08a283080dc27fff,
            0x08a283080ddb7fff, 0x08a283080dc07fff, 0x08a283080dd8ffff,
            0x08a283080dd5ffff, 0x08a283080dc4ffff, 0x08a283080dd47fff,
            0x08a283080dce7fff, 0x08a283080dd1ffff, 0x08a283080dceffff,
            0x08a283080dc6ffff, 0x08a283080dc87fff, 0x08a283080dcaffff,
            0x08a283080dd2ffff, 0x08a283080dcd7fff, 0x08a283080dd9ffff,
            0x08a283080dd6ffff, 0x08a283080dcc7fff, 0x08a283080dca7fff,
            0x08a283080dccffff, 0x08a283080dd77fff, 0x08a283080dc97fff,
            0x08a283080dd4ffff, 0x08a283080dd97fff, 0x08a283080dc37fff,
            0x08a283080dc8ffff, 0x08a283080dcb7fff, 0x08a283080dcf7fff,
            0x08a283080dd87fff, 0x08a283080dda7fff, 0x08a283080dc9ffff,
            0x08a283080dc77fff, 0x08a283080dc67fff, 0x08a283080dc57fff,
            0x08a283080ddaffff, 0x08a283080dd17fff, 0x08a283080dc17fff,
            0x08a283080dd57fff, 0x08a283080dc0ffff, 0x08a283080dd07fff,
            0x08a283080dc1ffff, 0x08a283080dd0ffff, 0x08a283080dc2ffff,
            0x08a283080dd67fff,
        ]
        @test Set(children) == Set(expected)
    end

    @testset "sameRes" begin
        h = H3X.H3Index(0x088283080ddfffff)
        (err, children) = H3X.cellToChildren(h, 8)
        @test err == H3X.E_SUCCESS
        @test children == [h]
    end

    @testset "childResTooCoarse" begin
        h = H3X.H3Index(0x088283080ddfffff)
        (err, _) = H3X.cellToChildren(h, 7)
        @test err == H3X.E_RES_DOMAIN
    end

    @testset "childResTooFine" begin
        h = H3X.H3Index(0x08f283080dcb0ae2)  # res 15 cell
        (err, _) = H3X.cellToChildren(h, H3X.MAX_H3_RES + 1)
        @test err == H3X.E_RES_DOMAIN
    end

    @testset "pentagonChildren" begin
        h = H3X.H3Index(0x081083ffffffffff)  # res 1 pentagon
        (err, children) = H3X.cellToChildren(h, 3)
        @test err == H3X.E_SUCCESS
        expected = H3X.H3Index[
            0x0830800fffffffff, 0x0830802fffffffff, 0x0830803fffffffff,
            0x0830804fffffffff, 0x0830805fffffffff, 0x0830806fffffffff,
            0x0830810fffffffff, 0x0830811fffffffff, 0x0830812fffffffff,
            0x0830813fffffffff, 0x0830814fffffffff, 0x0830815fffffffff,
            0x0830816fffffffff, 0x0830818fffffffff, 0x0830819fffffffff,
            0x083081afffffffff, 0x083081bfffffffff, 0x083081cfffffffff,
            0x083081dfffffffff, 0x083081efffffffff, 0x0830820fffffffff,
            0x0830821fffffffff, 0x0830822fffffffff, 0x0830823fffffffff,
            0x0830824fffffffff, 0x0830825fffffffff, 0x0830826fffffffff,
            0x0830828fffffffff, 0x0830829fffffffff, 0x083082afffffffff,
            0x083082bfffffffff, 0x083082cfffffffff, 0x083082dfffffffff,
            0x083082efffffffff, 0x0830830fffffffff, 0x0830831fffffffff,
            0x0830832fffffffff, 0x0830833fffffffff, 0x0830834fffffffff,
            0x0830835fffffffff, 0x0830836fffffffff,
        ]
        @test Set(children) == Set(expected)
    end
end

# ============================================================================
# testCellToChildrenSize
# ============================================================================
@testset "cellToChildrenSize" begin
    @testset "cellToChildrenSize_hexagon" begin
        h = H3X.H3Index(0x087283080dffffff)  # res 7 hexagon
        (err, sz) = H3X.cellToChildrenSize(h, 3)
        @test err == H3X.E_RES_DOMAIN
        (err2, sz2) = H3X.cellToChildrenSize(h, 7)
        @test err2 == H3X.E_SUCCESS
        @test sz2 == 1
        (err3, sz3) = H3X.cellToChildrenSize(h, 8)
        @test err3 == H3X.E_SUCCESS
        @test sz3 == 7
        (err4, sz4) = H3X.cellToChildrenSize(h, 9)
        @test err4 == H3X.E_SUCCESS
        @test sz4 == 7 * 7
    end

    @testset "cellToChildrenSize_pentagon" begin
        h = H3X.H3Index(0x0870800000ffffff)  # res 7 pentagon
        (err, _) = H3X.cellToChildrenSize(h, 3)
        @test err == H3X.E_RES_DOMAIN
        (err2, sz2) = H3X.cellToChildrenSize(h, 7)
        @test err2 == H3X.E_SUCCESS
        @test sz2 == 1
        (err3, sz3) = H3X.cellToChildrenSize(h, 8)
        @test err3 == H3X.E_SUCCESS
        @test sz3 == 6
        (err4, sz4) = H3X.cellToChildrenSize(h, 9)
        @test err4 == H3X.E_SUCCESS
        @test sz4 == (5 * 7) + (1 * 6)
    end

    @testset "cellToChildrenSize_largest_hexagon" begin
        h = H3X.H3Index(0x0806dfffffffffff)  # res 0 hexagon
        expected = Int64(4747561509943)  # 7^15
        (err, out) = H3X.cellToChildrenSize(h, 15)
        @test err == H3X.E_SUCCESS
        @test out == expected
    end

    @testset "cellToChildrenSize_largest_pentagon" begin
        h = H3X.H3Index(0x08009fffffffffff)  # res 0 pentagon
        expected = Int64(3956301258286)  # 1 + 5*(7^15 - 1)/6
        (err, out) = H3X.cellToChildrenSize(h, 15)
        @test err == H3X.E_SUCCESS
        @test out == expected
    end
end

# ============================================================================
# testCellToCenterChild
# ============================================================================
@testset "cellToCenterChild" begin
    baseHex = H3X.setH3Index(8, 4, 2)
    (_, baseCentroid) = H3X.cellToLatLng(baseHex)

    @testset "propertyTests" begin
        for res in 0:(H3X.MAX_H3_RES - 1)
            for childRes in (res + 1):H3X.MAX_H3_RES
                (err, h3Index) = H3X.latLngToCell(baseCentroid, res)
                @test err == H3X.E_SUCCESS
                (_, centroid) = H3X.cellToLatLng(h3Index)

                (err2, geoChild) = H3X.latLngToCell(centroid, childRes)
                @test err2 == H3X.E_SUCCESS
                (err3, centerChild) = H3X.cellToCenterChild(h3Index, childRes)
                @test err3 == H3X.E_SUCCESS

                @test centerChild == geoChild
                @test H3X.getResolution(centerChild) == childRes
                (err4, parent) = H3X.cellToParent(centerChild, res)
                @test err4 == H3X.E_SUCCESS
                @test parent == h3Index
            end
        end
    end

    @testset "sameRes" begin
        res = H3X.getResolution(baseHex)
        (err, child) = H3X.cellToCenterChild(baseHex, res)
        @test err == H3X.E_SUCCESS
        @test child == baseHex
    end

    @testset "invalidInputs" begin
        res = H3X.getResolution(baseHex)
        (err, _) = H3X.cellToCenterChild(baseHex, res - 1)
        @test err == H3X.E_RES_DOMAIN
        (err2, _) = H3X.cellToCenterChild(baseHex, -1)
        @test err2 == H3X.E_RES_DOMAIN
        (err3, _) = H3X.cellToCenterChild(baseHex, H3X.MAX_H3_RES + 1)
        @test err3 == H3X.E_RES_DOMAIN
    end
end

# ============================================================================
# testCompactCells
# ============================================================================
@testset "compactCells" begin
    @testset "res0children" begin
        parent = H3X.setH3Index(0, 0, 0)
        (err, arrSize) = H3X.cellToChildrenSize(parent, 1)
        @test err == H3X.E_SUCCESS
        (err2, children) = H3X.cellToChildren(parent, 1)
        @test err2 == H3X.E_SUCCESS
        (err3, compressed) = H3X.compactCells(children)
        @test err3 == H3X.E_SUCCESS
        non_zero = filter(!=(H3X.H3_NULL), compressed)
        @test length(non_zero) == 1
        @test non_zero[1] == parent
    end

    @testset "res0" begin
        hexCount = H3X.NUM_BASE_CELLS
        res0Hexes = H3X.H3Index[H3X.setH3Index(0, i, 0) for i in 0:(hexCount - 1)]
        (err, compressed) = H3X.compactCells(res0Hexes)
        @test err == H3X.E_SUCCESS
        @test Set(filter(!=(H3X.H3_NULL), compressed)) == Set(res0Hexes)
    end

    @testset "compactCells_duplicate" begin
        numHex = 10
        someHexagons = H3X.H3Index[H3X.setH3Index(5, 0, 2) for _ in 1:numHex]
        (err, _) = H3X.compactCells(someHexagons)
        @test err == H3X.E_DUPLICATE_INPUT
    end

    @testset "compactCells_duplicateMinimum" begin
        h3 = H3X.setH3Index(10, 0, 2)
        (err, arrSize) = H3X.cellToChildrenSize(h3, 11)
        @test err == H3X.E_SUCCESS
        (err2, children) = H3X.cellToChildren(h3, 11)
        @test err2 == H3X.E_SUCCESS
        push!(children, children[1])
        (err3, _) = H3X.compactCells(children)
        @test err3 == H3X.E_DUPLICATE_INPUT
    end

    @testset "compactCells_duplicatePentagonLimit" begin
        h3 = H3X.setH3Index(10, 4, 0)
        (err, arrSize) = H3X.cellToChildrenSize(h3, 11)
        @test err == H3X.E_SUCCESS
        (err2, children) = H3X.cellToChildren(h3, 11)
        @test err2 == H3X.E_SUCCESS
        (err3, centerChild) = H3X.cellToCenterChild(h3, 11)
        @test err3 == H3X.E_SUCCESS
        push!(children, centerChild)
        (err4, _) = H3X.compactCells(children)
        @test err4 == H3X.E_DUPLICATE_INPUT
    end

    @testset "compactCells_empty" begin
        (err, result) = H3X.compactCells(H3X.H3Index[])
        @test err == H3X.E_SUCCESS
    end

    @testset "compactCells_disparate" begin
        numHex = 7
        disparate = H3X.H3Index[H3X.setH3Index(1, i, H3X.CENTER_DIGIT) for i in 0:(numHex - 1)]
        (err, output) = H3X.compactCells(disparate)
        @test err == H3X.E_SUCCESS
        @test Set(output) == Set(disparate)
    end

    @testset "someHexagon" begin
        origin = H3X.setH3Index(1, 5, 0)
        (err, childrenSz) = H3X.uncompactCellsSize(H3X.H3Index[origin], 2)
        @test err == H3X.E_SUCCESS
        (err2, children) = H3X.uncompactCells(H3X.H3Index[origin], 2)
        @test err2 == H3X.E_SUCCESS

        (err3, result) = H3X.compactCells(children)
        @test err3 == H3X.E_SUCCESS
        non_zero = filter(!=(H3X.H3_NULL), result)
        @test length(non_zero) == 1
        @test non_zero[1] == origin
    end

    @testset "pentagon" begin
        pentagon = H3X.setH3Index(1, 4, 0)
        (err, childrenSz) = H3X.uncompactCellsSize(H3X.H3Index[pentagon], 2)
        @test err == H3X.E_SUCCESS
        (err2, children) = H3X.uncompactCells(H3X.H3Index[pentagon], 2)
        @test err2 == H3X.E_SUCCESS

        (err3, result) = H3X.compactCells(children)
        @test err3 == H3X.E_SUCCESS
        non_zero = filter(!=(H3X.H3_NULL), result)
        @test length(non_zero) == 1
        @test non_zero[1] == pentagon
    end

    @testset "large_uncompact_size_hexagon" begin
        cells = H3X.H3Index[H3X.H3Index(0x0806dfffffffffff)]  # res 0 hexagon
        expected = Int64(4747561509943)  # 7^15
        (err, out) = H3X.uncompactCellsSize(cells, 15)
        @test err == H3X.E_SUCCESS
        @test out == expected
    end

    @testset "large_uncompact_size_pentagon" begin
        cells = H3X.H3Index[H3X.H3Index(0x08009fffffffffff)]  # res 0 pentagon
        expected = Int64(3956301258286)  # 1 + 5*(7^15 - 1)/6
        (err, out) = H3X.uncompactCellsSize(cells, 15)
        @test err == H3X.E_SUCCESS
        @test out == expected
    end

    @testset "uncompactCells_wrongRes" begin
        someHexagons = H3X.H3Index[H3X.setH3Index(5, i, 0) for i in 0:2]

        (err, _) = H3X.uncompactCellsSize(someHexagons, 4)
        @test err == H3X.E_RES_MISMATCH
        (err2, _) = H3X.uncompactCellsSize(someHexagons, -1)
        @test err2 == H3X.E_RES_MISMATCH
        (err3, _) = H3X.uncompactCellsSize(someHexagons, H3X.MAX_H3_RES + 1)
        @test err3 == H3X.E_RES_MISMATCH
    end

    @testset "uncompactCells_empty" begin
        (err, sz) = H3X.uncompactCellsSize(H3X.H3Index[], 0)
        @test err == H3X.E_SUCCESS
        @test sz == 0
        (err2, result) = H3X.uncompactCells(H3X.H3Index[], 0)
        @test err2 == H3X.E_SUCCESS
    end

    @testset "uncompactCells_onlyZero" begin
        origin = H3X.H3Index[H3X.H3_NULL]
        (err, childrenSz) = H3X.uncompactCellsSize(origin, 2)
        @test err == H3X.E_SUCCESS
        (err2, children) = H3X.uncompactCells(origin, 2)
        @test err2 == H3X.E_SUCCESS
    end
end

# ============================================================================
# testPentagonIndexes
# ============================================================================
@testset "getPentagons" begin
    @testset "propertyTests" begin
        expectedCount = H3X.pentagonCount()
        for res in 0:15
            (err, h3Indexes) = H3X.getPentagons(res)
            @test err == H3X.E_SUCCESS

            numFound = 0
            for (i, h3Index) in enumerate(h3Indexes)
                if h3Index != H3X.H3_NULL
                    numFound += 1
                    @test H3X.isValidCell(h3Index)
                    @test H3X.isPentagon(h3Index)
                    @test H3X.getResolution(h3Index) == res

                    for j in (i + 1):length(h3Indexes)
                        @test h3Indexes[j] != h3Index
                    end
                end
            end
            @test numFound == expectedCount
        end
    end

    @testset "getPentagonsInvalid" begin
        (err, _) = H3X.getPentagons(16)
        @test err == H3X.E_RES_DOMAIN
        (err2, _) = H3X.getPentagons(100)
        @test err2 == H3X.E_RES_DOMAIN
        (err3, _) = H3X.getPentagons(-1)
        @test err3 == H3X.E_RES_DOMAIN
    end

    @testset "invalidPentagons" begin
        @test !H3X.isPentagon(H3X.H3Index(0))
        @test !H3X.isPentagon(H3X.H3Index(0x7fffffffffffffff))
    end
end

# ============================================================================
# testDescribeH3Error
# ============================================================================
@testset "describeH3Error" begin
    @testset "noError" begin
        @test H3X.describeH3Error(H3X.E_SUCCESS) == "Success"
    end

    @testset "invalidCell" begin
        @test H3X.describeH3Error(H3X.E_CELL_INVALID) == "H3Index cell argument was not valid"
    end

    @testset "allErrorCodes" begin
        for e in instances(H3X.H3Error)
            desc = H3X.describeH3Error(e)
            @test desc isa String
            @test !isempty(desc)
        end
    end
end

# ============================================================================
# testIndexDigits
# ============================================================================
@testset "indexDigits" begin
    @testset "getIndexDigitForCell" begin
        anywhere = H3X.LatLng(0.0, 0.0)
        local h
        for resCell in 0:H3X.MAX_H3_RES
            (err, h) = H3X.latLngToCell(anywhere, resCell)
            @test err == H3X.E_SUCCESS
            for resDigit in 1:H3X.MAX_H3_RES
                (err2, digit) = H3X.getIndexDigit(h, resDigit)
                @test err2 == H3X.E_SUCCESS
                if resDigit <= resCell
                    @test digit >= Int(H3X.CENTER_DIGIT) && digit < Int(H3X.INVALID_DIGIT)
                else
                    @test digit == Int(H3X.INVALID_DIGIT)
                end
            end
        end

        (err_neg, _) = H3X.getIndexDigit(h, -1)
        @test err_neg == H3X.E_RES_DOMAIN
        (err_zero, _) = H3X.getIndexDigit(h, 0)
        @test err_zero == H3X.E_RES_DOMAIN
        (err_high, _) = H3X.getIndexDigit(h, 16)
        @test err_high == H3X.E_RES_DOMAIN
    end

    @testset "getIndexDigitForSetCell" begin
        for expectedDigit in Int(H3X.CENTER_DIGIT):(Int(H3X.INVALID_DIGIT) - 1)
            for resCell in 0:H3X.MAX_H3_RES
                h = H3X.setH3Index(resCell, 0, expectedDigit)
                for resDigit in 1:H3X.MAX_H3_RES
                    (err, digit) = H3X.getIndexDigit(h, resDigit)
                    @test err == H3X.E_SUCCESS
                    if resDigit <= resCell
                        @test digit == expectedDigit
                    else
                        @test digit == Int(H3X.INVALID_DIGIT)
                    end
                end
            end
        end
    end
end

# ============================================================================
# testLatLngToCell (basic validation - the C version reads from stdin)
# ============================================================================
@testset "latLngToCell" begin
    @testset "roundtrip_basic" begin
        for res in 0:H3X.MAX_H3_RES
            g = H3X.LatLng(0.0, 0.0)
            (err, h1) = H3X.latLngToCell(g, res)
            @test err == H3X.E_SUCCESS
            (err2, g2) = H3X.cellToLatLng(h1)
            @test err2 == H3X.E_SUCCESS
            (err3, h2) = H3X.latLngToCell(g2, res)
            @test err3 == H3X.E_SUCCESS
            @test h1 == h2
        end
    end
end

# ============================================================================
# testCellToLatLng (basic validation - the C version reads from stdin)
# ============================================================================
@testset "cellToLatLng" begin
    @testset "roundtrip_h_to_latlng_to_h" begin
        for res in 0:H3X.MAX_H3_RES
            g1 = H3X.LatLng(0.0, 0.0)
            (err, h1) = H3X.latLngToCell(g1, res)
            @test err == H3X.E_SUCCESS
            (err2, g2) = H3X.cellToLatLng(h1)
            @test err2 == H3X.E_SUCCESS
            (err3, h2) = H3X.latLngToCell(g2, res)
            @test err3 == H3X.E_SUCCESS
            @test h1 == h2
        end
    end
end

# ============================================================================
# testCellToBoundary (basic validation - the C version reads from stdin)
# ============================================================================
@testset "cellToBoundary" begin
    @testset "basic_hex_boundary" begin
        g = H3X.LatLng(0.0, 0.0)
        for res in 0:5
            (err, h) = H3X.latLngToCell(g, res)
            @test err == H3X.E_SUCCESS
            (err2, b) = H3X.cellToBoundary(h)
            @test err2 == H3X.E_SUCCESS
            if H3X.isPentagon(h)
                @test b.numVerts == 5
            else
                @test b.numVerts >= 6 && b.numVerts <= 7
            end
        end
    end
end
