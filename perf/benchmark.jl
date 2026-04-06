#=
Benchmark gridPathCells across a wide range of inputs.

Compares three implementations:
  - cube:   FastH3.gridPathCells (cube-coordinate interpolation)
  - gc:     FastH3.FastH3Extension.gridPathCells (great-circle walk)
  - hybrid: FastH3.FastH3Extension.gridPathCellsHybrid (piecewise cube + GC fallback)

Usage:
  julia --project=perf perf/benchmark.jl
  julia --project=perf perf/benchmark.jl --resolution 5
  julia --project=perf perf/benchmark.jl --distance 10
  julia --project=perf perf/benchmark.jl --scenario pent
=#

include("common.jl")
import BenchmarkTools
import Printf

function _fmt_time(ns::Float64)
    ns < 1e3 && return Printf.@sprintf("%.1f ns", ns)
    ns < 1e6 && return Printf.@sprintf("%.1f μs", ns / 1e3)
    ns < 1e9 && return Printf.@sprintf("%.1f ms", ns / 1e6)
    return Printf.@sprintf("%.2f s", ns / 1e9)
end

function main()
    args = parse_perf_args(; description="Benchmark gridPathCells (cube vs gc vs hybrid)")
    scenarios = filter_scenarios(build_scenarios(), args)

    if isempty(scenarios)
        println("No scenarios matched the filters.")
        return
    end

    println("gridPathCells benchmark — cube vs gc vs hybrid ($(length(scenarios)) scenarios)\n")
    hdr = string(
        rpad("scenario", 22), " │ ",
        rpad("res", 4), "  ",
        rpad("dist (hops)", 11), " │ ",
        rpad("cells(cube)", 11), " ",
        rpad("cells(gc)", 11), " ",
        rpad("cells(hyb)", 11), " │ ",
        rpad("min(cube)", 12), " ",
        rpad("min(gc)", 12), " ",
        rpad("min(hyb)", 12), " │ ",
        rpad("alloc(cube)", 11), " ",
        rpad("alloc(gc)", 11), " ",
        rpad("alloc(hyb)", 11), " │ ",
        rpad("mem(cube)", 10), " ",
        rpad("mem(gc)", 10), " ",
        rpad("mem(hyb)", 10), " │",
    )
    println(hdr)
    println(replace(hdr, r"[^│\n]" => "─"))

    Ext = FastH3.FastH3Extension

    for sc in scenarios
        FastH3.gridPathCells(sc.start_cell, sc.end_cell)
        Ext.gridPathCells(sc.start_cell, sc.end_cell)
        Ext.gridPathCellsHybrid(sc.start_cell, sc.end_cell)

        b_cube = BenchmarkTools.@benchmark FastH3.gridPathCells($(sc.start_cell), $(sc.end_cell))
        err_c, path_c = FastH3.gridPathCells(sc.start_cell, sc.end_cell)
        n_cube = err_c == FastH3.E_SUCCESS ? length(path_c) : 0

        b_gc = BenchmarkTools.@benchmark $(Ext).gridPathCells($(sc.start_cell), $(sc.end_cell))
        path_gc = Ext.gridPathCells(sc.start_cell, sc.end_cell)
        n_gc = length(path_gc)

        b_hyb = BenchmarkTools.@benchmark $(Ext).gridPathCellsHybrid($(sc.start_cell), $(sc.end_cell))
        path_hyb = Ext.gridPathCellsHybrid(sc.start_cell, sc.end_cell)
        n_hyb = length(path_hyb)

        println(
            rpad(sc.name, 22), " │ ",
            rpad(string(sc.resolution), 4), "  ",
            rpad(string(sc.distance), 11), " │ ",
            rpad(string(n_cube), 11), " ",
            rpad(string(n_gc), 11), " ",
            rpad(string(n_hyb), 11), " │ ",
            rpad(_fmt_time(minimum(b_cube.times)), 12), " ",
            rpad(_fmt_time(minimum(b_gc.times)), 12), " ",
            rpad(_fmt_time(minimum(b_hyb.times)), 12), " │ ",
            rpad(string(b_cube.allocs), 11), " ",
            rpad(string(b_gc.allocs), 11), " ",
            rpad(string(b_hyb.allocs), 11), " │ ",
            rpad(BenchmarkTools.prettymemory(b_cube.memory), 10), " ",
            rpad(BenchmarkTools.prettymemory(b_gc.memory), 10), " ",
            rpad(BenchmarkTools.prettymemory(b_hyb.memory), 10), " │",
        )
    end
end

main()
