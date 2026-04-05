# Tests are included from runtests.jl which loads FastH3

function testDecreasingFunction(f)
    last = 0.0
    for i in FastH3.MAX_H3_RES:-1:0
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
        degs = FastH3.radsToDegs(originalRads)
        rads = FastH3.degsToRads(degs)
        @test abs(rads - originalRads) < FastH3.EPSILON_RAD
    end

    @testset "distanceRads" begin
        p1 = FastH3.setGeoDegs(10.0, 10.0)
        p2 = FastH3.setGeoDegs(0.0, 10.0)
        @test FastH3.greatCircleDistanceRads(p1, p1) < FastH3.EPSILON_RAD * 1000
        @test abs(FastH3.greatCircleDistanceRads(p1, p2) - FastH3.degsToRads(10.0)) < FastH3.EPSILON_RAD * 1000
    end

    @testset "distanceRads_wrappedLongitude" begin
        negativeLongitude = LatLng(0.0, -(FastH3.M_PI + FastH3.M_PI_2))
        zero = LatLng(0.0, 0.0)
        @test abs(FastH3.M_PI_2 - FastH3.greatCircleDistanceRads(negativeLongitude, zero)) < FastH3.EPSILON_RAD
        @test abs(FastH3.M_PI_2 - FastH3.greatCircleDistanceRads(zero, negativeLongitude)) < FastH3.EPSILON_RAD
    end

    @testset "doubleConstants" begin
        testDecreasingFunction(FastH3.getHexagonAreaAvgKm2)
        testDecreasingFunction(FastH3.getHexagonAreaAvgM2)
        testDecreasingFunction(FastH3.getHexagonEdgeLengthAvgKm)
        testDecreasingFunction(FastH3.getHexagonEdgeLengthAvgM)
    end

    @testset "doubleConstantsErrors" begin
        @test first(FastH3.getHexagonAreaAvgKm2(-1)) == E_RES_DOMAIN
        @test first(FastH3.getHexagonAreaAvgKm2(16)) == E_RES_DOMAIN
        @test first(FastH3.getHexagonAreaAvgM2(-1)) == E_RES_DOMAIN
        @test first(FastH3.getHexagonAreaAvgM2(16)) == E_RES_DOMAIN
        @test first(FastH3.getHexagonEdgeLengthAvgKm(-1)) == E_RES_DOMAIN
        @test first(FastH3.getHexagonEdgeLengthAvgKm(16)) == E_RES_DOMAIN
        @test first(FastH3.getHexagonEdgeLengthAvgM(-1)) == E_RES_DOMAIN
        @test first(FastH3.getHexagonEdgeLengthAvgM(16)) == E_RES_DOMAIN
    end

    @testset "intConstants" begin
        last = Int64(0)
        for i in 0:FastH3.MAX_H3_RES
            (err, next) = FastH3.getNumCells(i)
            @test err == E_SUCCESS
            @test next > last
            last = next
        end
    end

    @testset "intConstantsErrors" begin
        @test first(FastH3.getNumCells(-1)) == E_RES_DOMAIN
        @test first(FastH3.getNumCells(16)) == E_RES_DOMAIN
    end

    @testset "numHexagons" begin
        expected = Int64[
            122, 842, 5882, 41162, 288122, 2016842,
            14117882, 98825162, 691776122, 4842432842,
            33897029882, 237279209162, 1660954464122,
            11626681248842, 81386768741882, 569707381193162,
        ]
        for r in 0:FastH3.MAX_H3_RES
            (err, num) = FastH3.getNumCells(r)
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
        @test FastH3.geoAlmostEqualThreshold(a, b, eps(Float64))

        b = LatLng(15.00001, 10.00002)
        @test FastH3.geoAlmostEqualThreshold(a, b, 0.0001)

        b = LatLng(15.00001, 10.0)
        @test !FastH3.geoAlmostEqualThreshold(a, b, 0.000001)

        b = LatLng(15.0, 10.00001)
        @test !FastH3.geoAlmostEqualThreshold(a, b, 0.000001)
    end

    @testset "constrainLatLng" begin
        @test FastH3.constrainLat(0.0) == 0.0
        @test FastH3.constrainLat(1.0) == 1.0
        @test FastH3.constrainLat(FastH3.M_PI_2) == FastH3.M_PI_2
        @test FastH3.constrainLat(FastH3.M_PI) == 0.0
        @test FastH3.constrainLat(FastH3.M_PI + 1.0) == 1.0
        @test FastH3.constrainLat(2.0 * FastH3.M_PI + 1.0) == 1.0

        @test FastH3.constrainLng(0.0) == 0.0
        @test FastH3.constrainLng(1.0) == 1.0
        @test FastH3.constrainLng(FastH3.M_PI) == FastH3.M_PI
        @test FastH3.constrainLng(2.0 * FastH3.M_PI) == 0.0
        @test FastH3.constrainLng(3.0 * FastH3.M_PI) == FastH3.M_PI
        @test FastH3.constrainLng(4.0 * FastH3.M_PI) == 0.0
    end

    @testset "_geoAzDistanceRads_noop" begin
        start = LatLng(15.0, 10.0)
        expected = LatLng(15.0, 10.0)
        out = FastH3._geoAzDistanceRads(start, 0.0, 0.0)
        @test FastH3.geoAlmostEqual(expected, out)
    end

    @testset "_geoAzDistanceRads_dueNorthSouth" begin
        start = FastH3.setGeoDegs(45.0, 1.0)
        expected = FastH3.setGeoDegs(90.0, 0.0)
        out = FastH3._geoAzDistanceRads(start, 0.0, FastH3.degsToRads(45.0))
        @test FastH3.geoAlmostEqual(expected, out)

        start = FastH3.setGeoDegs(45.0, 1.0)
        expected = FastH3.setGeoDegs(270.0, 1.0)
        out = FastH3._geoAzDistanceRads(start, 0.0, FastH3.degsToRads(225.0))
        @test FastH3.geoAlmostEqual(expected, out)

        start = FastH3.setGeoDegs(-45.0, 2.0)
        expected = FastH3.setGeoDegs(-90.0, 0.0)
        out = FastH3._geoAzDistanceRads(start, FastH3.degsToRads(180.0), FastH3.degsToRads(45.0))
        @test FastH3.geoAlmostEqual(expected, out)

        start = FastH3.setGeoDegs(-45.0, 10.0)
        expected = FastH3.setGeoDegs(-10.0, 10.0)
        out = FastH3._geoAzDistanceRads(start, 0.0, FastH3.degsToRads(35.0))
        @test FastH3.geoAlmostEqual(expected, out)
    end

    @testset "_geoAzDistanceRads_poleToPole" begin
        start = FastH3.setGeoDegs(90.0, 0.0)
        expected = FastH3.setGeoDegs(-90.0, 0.0)
        out = FastH3._geoAzDistanceRads(start, FastH3.degsToRads(12.0), FastH3.degsToRads(180.0))
        @test FastH3.geoAlmostEqual(expected, out)

        start = FastH3.setGeoDegs(-90.0, 0.0)
        expected = FastH3.setGeoDegs(90.0, 0.0)
        out = FastH3._geoAzDistanceRads(start, FastH3.degsToRads(34.0), FastH3.degsToRads(180.0))
        @test FastH3.geoAlmostEqual(expected, out)
    end

    @testset "_geoAzDistanceRads_invertible" begin
        start = FastH3.setGeoDegs(15.0, 10.0)
        azimuth = FastH3.degsToRads(20.0)
        degrees180 = FastH3.degsToRads(180.0)
        distance = FastH3.degsToRads(15.0)

        out = FastH3._geoAzDistanceRads(start, azimuth, distance)
        @test abs(FastH3.greatCircleDistanceRads(start, out) - distance) < FastH3.EPSILON_RAD

        start2 = out
        out2 = FastH3._geoAzDistanceRads(start2, azimuth + degrees180, distance)
        @test FastH3.greatCircleDistanceRads(start, out2) < 0.01
    end
end

# ---------------------------------------------------------------------------- #
#  testMathExtensionsInternal.c
# ---------------------------------------------------------------------------- #

@testset "mathExtensionsInternal" begin
    @testset "_ipow" begin
        @test FastH3._ipow(7, 0) == 1
        @test FastH3._ipow(7, 1) == 7
        @test FastH3._ipow(7, 2) == 49
        @test FastH3._ipow(1, 20) == 1
        @test FastH3._ipow(2, 5) == 32
    end

    @testset "subOverflows" begin
        @test !FastH3.SUB_INT32S_OVERFLOWS(Int32(0), Int32(0))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32), Int32(0))
        @test  FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32), Int32(1))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32), Int32(-1))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(0))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(1))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-1))
        @test  FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(2))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-2))
        @test !FastH3.SUB_INT32S_OVERFLOWS(Int32(100), Int32(10))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32), Int32(0))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32), Int32(1))
        @test  FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32), Int32(-1))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(1))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-1))
        @test  FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-2))
        @test  FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-2))
        @test  FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32), typemax(Int32))
        @test  FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32), typemin(Int32))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemin(Int32), typemin(Int32))
        @test !FastH3.SUB_INT32S_OVERFLOWS(typemax(Int32), typemax(Int32))
        @test !FastH3.SUB_INT32S_OVERFLOWS(Int32(-1), Int32(0))
        @test !FastH3.SUB_INT32S_OVERFLOWS(Int32(-1), Int32(10))
        @test !FastH3.SUB_INT32S_OVERFLOWS(Int32(-1), Int32(-10))
        @test !FastH3.SUB_INT32S_OVERFLOWS(Int32(-1), typemax(Int32))
        @test  FastH3.SUB_INT32S_OVERFLOWS(Int32(-2), typemax(Int32))
        @test !FastH3.SUB_INT32S_OVERFLOWS(Int32(-1), typemin(Int32))
        @test  FastH3.SUB_INT32S_OVERFLOWS(Int32(0), typemin(Int32))
    end

    @testset "addOverflows" begin
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(0), Int32(0))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32), Int32(0))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32), Int32(1))
        @test  FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32), Int32(-1))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(0))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(1))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-1))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(2))
        @test  FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32) + Int32(1), Int32(-2))
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(100), Int32(10))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32), Int32(0))
        @test  FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32), Int32(1))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32), Int32(-1))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(1))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-1))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(-2))
        @test  FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32) - Int32(1), Int32(2))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32), typemax(Int32))
        @test !FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32), typemin(Int32))
        @test  FastH3.ADD_INT32S_OVERFLOWS(typemax(Int32), typemax(Int32))
        @test  FastH3.ADD_INT32S_OVERFLOWS(typemin(Int32), typemin(Int32))
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(-1), Int32(0))
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(-1), Int32(10))
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(-1), Int32(-10))
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(-1), typemax(Int32))
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(-2), typemax(Int32))
        @test  FastH3.ADD_INT32S_OVERFLOWS(Int32(-1), typemin(Int32))
        @test !FastH3.ADD_INT32S_OVERFLOWS(Int32(0), typemin(Int32))
    end
end
