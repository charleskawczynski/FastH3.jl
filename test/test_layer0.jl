# Tests are included from runtests.jl which loads H3X

function testDecreasingFunction(f)
    last = 0.0
    for i in H3X.MAX_H3_RES:-1:0
        (err, next) = f(i)
        @test err == E_SUCCESS
        @test next > last
        last = next
    end
end

# ---------------------------------------------------------------------------- #
#  testLatLng.c
# ---------------------------------------------------------------------------- #

@testset "latLng" begin
    @testset "radsToDegs" begin
        originalRads = 1.0
        degs = H3X.radsToDegs(originalRads)
        rads = H3X.degsToRads(degs)
        @test abs(rads - originalRads) < H3X.EPSILON_RAD
    end

    @testset "distanceRads" begin
        p1 = H3X.setGeoDegs(10.0, 10.0)
        p2 = H3X.setGeoDegs(0.0, 10.0)
        @test H3X.greatCircleDistanceRads(p1, p1) < H3X.EPSILON_RAD * 1000
        @test abs(H3X.greatCircleDistanceRads(p1, p2) - H3X.degsToRads(10.0)) < H3X.EPSILON_RAD * 1000
    end

    @testset "distanceRads_wrappedLongitude" begin
        negativeLongitude = LatLng(0.0, -(H3X.M_PI + H3X.M_PI_2))
        zero = LatLng(0.0, 0.0)
        @test abs(H3X.M_PI_2 - H3X.greatCircleDistanceRads(negativeLongitude, zero)) < H3X.EPSILON_RAD
        @test abs(H3X.M_PI_2 - H3X.greatCircleDistanceRads(zero, negativeLongitude)) < H3X.EPSILON_RAD
    end

    @testset "doubleConstants" begin
        testDecreasingFunction(H3X.getHexagonAreaAvgKm2)
        testDecreasingFunction(H3X.getHexagonAreaAvgM2)
        testDecreasingFunction(H3X.getHexagonEdgeLengthAvgKm)
        testDecreasingFunction(H3X.getHexagonEdgeLengthAvgM)
    end

    @testset "doubleConstantsErrors" begin
        @test first(H3X.getHexagonAreaAvgKm2(-1)) == E_RES_DOMAIN
        @test first(H3X.getHexagonAreaAvgKm2(16)) == E_RES_DOMAIN
        @test first(H3X.getHexagonAreaAvgM2(-1)) == E_RES_DOMAIN
        @test first(H3X.getHexagonAreaAvgM2(16)) == E_RES_DOMAIN
        @test first(H3X.getHexagonEdgeLengthAvgKm(-1)) == E_RES_DOMAIN
        @test first(H3X.getHexagonEdgeLengthAvgKm(16)) == E_RES_DOMAIN
        @test first(H3X.getHexagonEdgeLengthAvgM(-1)) == E_RES_DOMAIN
        @test first(H3X.getHexagonEdgeLengthAvgM(16)) == E_RES_DOMAIN
    end

    @testset "intConstants" begin
        last = Int64(0)
        for i in 0:H3X.MAX_H3_RES
            (err, next) = H3X.getNumCells(i)
            @test err == E_SUCCESS
            @test next > last
            last = next
        end
    end

    @testset "intConstantsErrors" begin
        @test first(H3X.getNumCells(-1)) == E_RES_DOMAIN
        @test first(H3X.getNumCells(16)) == E_RES_DOMAIN
    end

    @testset "numHexagons" begin
        expected = Int64[
            122, 842, 5882, 41162, 288122, 2016842,
            14117882, 98825162, 691776122, 4842432842,
            33897029882, 237279209162, 1660954464122,
            11626681248842, 81386768741882, 569707381193162,
        ]
        for r in 0:H3X.MAX_H3_RES
            (err, num) = H3X.getNumCells(r)
            @test err == E_SUCCESS
            @test num == expected[r + 1]
        end
    end
end

# ---------------------------------------------------------------------------- #
#  testLatLngInternal.c
# ---------------------------------------------------------------------------- #

@testset "latLngInternal" begin
    @testset "geoAlmostEqualThreshold" begin
        a = LatLng(15.0, 10.0)
        b = LatLng(15.0, 10.0)
        @test H3X.geoAlmostEqualThreshold(a, b, eps(Float64))

        b = LatLng(15.00001, 10.00002)
        @test H3X.geoAlmostEqualThreshold(a, b, 0.0001)

        b = LatLng(15.00001, 10.0)
        @test !H3X.geoAlmostEqualThreshold(a, b, 0.000001)

        b = LatLng(15.0, 10.00001)
        @test !H3X.geoAlmostEqualThreshold(a, b, 0.000001)
    end

    @testset "constrainLatLng" begin
        @test H3X.constrainLat(0.0) == 0.0
        @test H3X.constrainLat(1.0) == 1.0
        @test H3X.constrainLat(H3X.M_PI_2) == H3X.M_PI_2
        @test H3X.constrainLat(H3X.M_PI) == 0.0
        @test H3X.constrainLat(H3X.M_PI + 1.0) == 1.0
        @test H3X.constrainLat(2.0 * H3X.M_PI + 1.0) == 1.0

        @test H3X.constrainLng(0.0) == 0.0
        @test H3X.constrainLng(1.0) == 1.0
        @test H3X.constrainLng(H3X.M_PI) == H3X.M_PI
        @test H3X.constrainLng(2.0 * H3X.M_PI) == 0.0
        @test H3X.constrainLng(3.0 * H3X.M_PI) == H3X.M_PI
        @test H3X.constrainLng(4.0 * H3X.M_PI) == 0.0
    end

    @testset "_geoAzDistanceRads_noop" begin
        start = LatLng(15.0, 10.0)
        expected = LatLng(15.0, 10.0)
        out = H3X._geoAzDistanceRads(start, 0.0, 0.0)
        @test H3X.geoAlmostEqual(expected, out)
    end

    @testset "_geoAzDistanceRads_dueNorthSouth" begin
        start = H3X.setGeoDegs(45.0, 1.0)
        expected = H3X.setGeoDegs(90.0, 0.0)
        out = H3X._geoAzDistanceRads(start, 0.0, H3X.degsToRads(45.0))
        @test H3X.geoAlmostEqual(expected, out)

        start = H3X.setGeoDegs(45.0, 1.0)
        expected = H3X.setGeoDegs(270.0, 1.0)
        out = H3X._geoAzDistanceRads(start, 0.0, H3X.degsToRads(225.0))
        @test H3X.geoAlmostEqual(expected, out)

        start = H3X.setGeoDegs(-45.0, 2.0)
        expected = H3X.setGeoDegs(-90.0, 0.0)
        out = H3X._geoAzDistanceRads(start, H3X.degsToRads(180.0), H3X.degsToRads(45.0))
        @test H3X.geoAlmostEqual(expected, out)

        start = H3X.setGeoDegs(-45.0, 10.0)
        expected = H3X.setGeoDegs(-10.0, 10.0)
        out = H3X._geoAzDistanceRads(start, 0.0, H3X.degsToRads(35.0))
        @test H3X.geoAlmostEqual(expected, out)
    end

    @testset "_geoAzDistanceRads_poleToPole" begin
        start = H3X.setGeoDegs(90.0, 0.0)
        expected = H3X.setGeoDegs(-90.0, 0.0)
        out = H3X._geoAzDistanceRads(start, H3X.degsToRads(12.0), H3X.degsToRads(180.0))
        @test H3X.geoAlmostEqual(expected, out)

        start = H3X.setGeoDegs(-90.0, 0.0)
        expected = H3X.setGeoDegs(90.0, 0.0)
        out = H3X._geoAzDistanceRads(start, H3X.degsToRads(34.0), H3X.degsToRads(180.0))
        @test H3X.geoAlmostEqual(expected, out)
    end

    @testset "_geoAzDistanceRads_invertible" begin
        start = H3X.setGeoDegs(15.0, 10.0)
        azimuth = H3X.degsToRads(20.0)
        degrees180 = H3X.degsToRads(180.0)
        distance = H3X.degsToRads(15.0)

        out = H3X._geoAzDistanceRads(start, azimuth, distance)
        @test abs(H3X.greatCircleDistanceRads(start, out) - distance) < H3X.EPSILON_RAD

        start2 = out
        out2 = H3X._geoAzDistanceRads(start2, azimuth + degrees180, distance)
        @test H3X.greatCircleDistanceRads(start, out2) < 0.01
    end
end

# ---------------------------------------------------------------------------- #
#  testMathExtensionsInternal.c
# ---------------------------------------------------------------------------- #

@testset "mathExtensionsInternal" begin
    @testset "_ipow" begin
        @test H3X._ipow(7, 0) == 1
        @test H3X._ipow(7, 1) == 7
        @test H3X._ipow(7, 2) == 49
        @test H3X._ipow(1, 20) == 1
        @test H3X._ipow(2, 5) == 32
    end

    @testset "subOverflows" begin
        @test !H3X.SUB_INT32S_OVERFLOWS(Int32(0), Int32(0))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemin(Int32), Int32(0))
        @test  H3X.SUB_INT32S_OVERFLOWS(typemin(Int32), Int32(1))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemin(Int32), Int32(-1))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(0))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(1))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-1))
        @test  H3X.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(2))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-2))
        @test !H3X.SUB_INT32S_OVERFLOWS(Int32(100), Int32(10))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemax(Int32), Int32(0))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemax(Int32), Int32(1))
        @test  H3X.SUB_INT32S_OVERFLOWS(typemax(Int32), Int32(-1))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(1))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-1))
        @test  H3X.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-2))
        @test  H3X.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-2))
        @test  H3X.SUB_INT32S_OVERFLOWS(typemin(Int32), typemax(Int32))
        @test  H3X.SUB_INT32S_OVERFLOWS(typemax(Int32), typemin(Int32))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemin(Int32), typemin(Int32))
        @test !H3X.SUB_INT32S_OVERFLOWS(typemax(Int32), typemax(Int32))
        @test !H3X.SUB_INT32S_OVERFLOWS(Int32(-1), Int32(0))
        @test !H3X.SUB_INT32S_OVERFLOWS(Int32(-1), Int32(10))
        @test !H3X.SUB_INT32S_OVERFLOWS(Int32(-1), Int32(-10))
        @test !H3X.SUB_INT32S_OVERFLOWS(Int32(-1), typemax(Int32))
        @test  H3X.SUB_INT32S_OVERFLOWS(Int32(-2), typemax(Int32))
        @test !H3X.SUB_INT32S_OVERFLOWS(Int32(-1), typemin(Int32))
        @test  H3X.SUB_INT32S_OVERFLOWS(Int32(0), typemin(Int32))
    end

    @testset "addOverflows" begin
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(0), Int32(0))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemin(Int32), Int32(0))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemin(Int32), Int32(1))
        @test  H3X.ADD_INT32S_OVERFLOWS(typemin(Int32), Int32(-1))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(0))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(1))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-1))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(2))
        @test  H3X.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-2))
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(100), Int32(10))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemax(Int32), Int32(0))
        @test  H3X.ADD_INT32S_OVERFLOWS(typemax(Int32), Int32(1))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemax(Int32), Int32(-1))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(1))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-1))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-2))
        @test  H3X.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(2))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemin(Int32), typemax(Int32))
        @test !H3X.ADD_INT32S_OVERFLOWS(typemax(Int32), typemin(Int32))
        @test  H3X.ADD_INT32S_OVERFLOWS(typemax(Int32), typemax(Int32))
        @test  H3X.ADD_INT32S_OVERFLOWS(typemin(Int32), typemin(Int32))
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(-1), Int32(0))
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(-1), Int32(10))
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(-1), Int32(-10))
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(-1), typemax(Int32))
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(-2), typemax(Int32))
        @test  H3X.ADD_INT32S_OVERFLOWS(Int32(-1), typemin(Int32))
        @test !H3X.ADD_INT32S_OVERFLOWS(Int32(0), typemin(Int32))
    end
end
