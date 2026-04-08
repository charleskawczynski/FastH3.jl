# Tests for `FastH3.FastH3Extension` (great-circle path, robust path, gridDistanceRobust).

const Ext = FastH3.FastH3Extension

@testset "FastH3Extension gridDistanceRobust" begin
    bc1 = FastH3.setH3Index(0, 15, 0)
    pent1 = FastH3.setH3Index(0, 4, 0)
    bc3 = FastH3.setH3Index(0, 31, 0)

    @testset "matches gridDistance for valid cells where core succeeds" begin
        err, d = FastH3.gridDistance(bc1, pent1)
        @test err == FastH3.E_SUCCESS
        err_r, d_r = Ext.gridDistanceRobust(bc1, pent1)
        @test err_r == FastH3.E_SUCCESS
        @test d_r == d
    end

    @testset "E_DOMAIN from gridDistance still yields robust distance" begin
        err, _ = FastH3.gridDistance(pent1, bc3)
        @test err == FastH3.E_DOMAIN
        err_r, d_r = Ext.gridDistanceRobust(pent1, bc3)
        @test err_r == FastH3.E_SUCCESS
        @test d_r == length(Ext.gridPathCellsRobust(pent1, bc3)) - 1
        @test d_r > 0
    end

    @testset "resolution and invalid indices match gridDistance errors" begin
        err_r, _ = Ext.gridDistanceRobust(
            FastH3.H3Index(0x0832830fffffffff),
            FastH3.H3Index(0x0822837fffffffff),
        )
        @test err_r == FastH3.E_RES_MISMATCH

        invalid = FastH3.H3Index(0xffffffffffffffff)
        err_r2, _ = Ext.gridDistanceRobust(invalid, invalid)
        @test err_r2 == FastH3.E_CELL_INVALID

        err_r3, _ = Ext.gridDistanceRobust(bc1, invalid)
        @test err_r3 == FastH3.E_RES_MISMATCH
    end

    @testset "directed edge indices defer to gridDistance" begin
        origin = FastH3.H3Index(0x0832830fffffffff)
        dest = FastH3.H3Index(0x0832834fffffffff)
        err_e, edge = FastH3.cellsToDirectedEdge(origin, dest)
        @test err_e == FastH3.E_SUCCESS
        err, d = FastH3.gridDistance(edge, origin)
        @test err == FastH3.E_SUCCESS
        err_r, d_r = Ext.gridDistanceRobust(edge, origin)
        @test err_r == FastH3.E_SUCCESS
        @test d_r == d
    end
end

@testset "FastH3Extension gridPathCellsRobust" begin
    @testset "cross-face endpoints and gridDistanceRobust" begin
        start_ = FastH3.H3Index(0x085285aa7fffffff)
        end_ = FastH3.H3Index(0x0851d9b1bfffffff)
        path = Ext.gridPathCellsRobust(start_, end_)
        @test !isempty(path)
        @test path[1] == start_
        @test path[end] == end_
        err_r, d_r = Ext.gridDistanceRobust(start_, end_)
        @test err_r == FastH3.E_SUCCESS
        @test d_r == length(path) - 1
    end

    @testset "same-face path matches cube and is a neighbor chain" begin
        start_ = FastH3.H3Index(0x0820807fffffffff)
        end_ = FastH3.H3Index(0x08208e7fffffffff)
        err_c, cube_path = FastH3.gridPathCells(start_, end_)
        @test err_c == FastH3.E_SUCCESS
        path = Ext.gridPathCellsRobust(start_, end_)
        @test path == cube_path
        for i in 2:length(path)
            err_n, is_n = FastH3.areNeighborCells(path[i], path[i - 1])
            @test err_n == FastH3.E_SUCCESS
            @test is_n
        end
    end
end

@testset "FastH3Extension gridPathCells (great-circle walk)" begin
    start_ = FastH3.H3Index(0x0820807fffffffff)
    end_ = FastH3.H3Index(0x08208e7fffffffff)
    path = Ext.gridPathCells(start_, end_)
    @test !isempty(path)
    @test path[1] == start_
    @test path[end] == end_
end
