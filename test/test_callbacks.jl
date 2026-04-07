@testset "Callback enumeration APIs" begin
    cell = H3Index(0x088283080ddfffff)
    neighbor = H3Index(0x088283080dc3ffff)
    cells = H3Index[cell]

    @testset "cellToChildren" begin
        err, vec = cellToChildren(cell, 8)
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = cellToChildren(cell, 8) do c
            push!(acc, c)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "uncompactCells" begin
        err, vec = uncompactCells(cells, 8)
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = uncompactCells(cells, 8) do c
            push!(acc, c)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "getPentagons" begin
        err, vec = getPentagons(5)
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = getPentagons(5) do h
            push!(acc, h)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "getRes0Cells" begin
        err, vec = getRes0Cells()
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = getRes0Cells() do h
            push!(acc, h)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "getIcosahedronFaces" begin
        err, vec = getIcosahedronFaces(cell)
        @assert err == E_SUCCESS
        acc = Int[]
        err2 = getIcosahedronFaces(cell) do f
            push!(acc, f)
        end
        @test err2 == E_SUCCESS
        want = filter(!=(-1), vec)
        @test acc == want
    end

    @testset "gridDiskDistancesUnsafe" begin
        err, v, d = gridDiskDistancesUnsafe(cell, 2)
        @assert err == E_SUCCESS
        accv = H3Index[]
        accd = Int[]
        err2 = gridDiskDistancesUnsafe(cell, 2) do c, r
            push!(accv, c)
            push!(accd, r)
        end
        @test err2 == E_SUCCESS
        @test accv == v
        @test accd == d
    end

    @testset "gridDiskUnsafe" begin
        err, vec = gridDiskUnsafe(cell, 2)
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = gridDiskUnsafe(cell, 2) do c
            push!(acc, c)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "gridRingUnsafe" begin
        err, vec = gridRingUnsafe(cell, 2)
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = gridRingUnsafe(cell, 2) do c
            push!(acc, c)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "gridPathCells" begin
        err, vec = gridPathCells(cell, neighbor)
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = gridPathCells(cell, neighbor) do c
            push!(acc, c)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "cellToVertexes" begin
        err, vec = cellToVertexes(cell)
        @assert err == E_SUCCESS
        acc = H3Index[]
        err2 = cellToVertexes(cell) do v
            push!(acc, v)
        end
        @test err2 == E_SUCCESS
        @test acc == vec
    end

    @testset "originToDirectedEdges" begin
        err_e, vec = originToDirectedEdges(cell)
        @assert err_e == E_SUCCESS
        acc = H3Index[]
        err3 = originToDirectedEdges(cell) do e
            push!(acc, e)
        end
        @test err3 == E_SUCCESS
        @test acc == vec
    end
end
