#=
Profile gridPathCells and write interactive flame graphs to perf/output/.

Supports profiling the cube-interpolation implementation, the great-circle
walk implementation, the hybrid implementation, or all (default).

Usage:
  julia --project=perf perf/flame.jl
  julia --project=perf perf/flame.jl --method gc
  julia --project=perf perf/flame.jl --method hybrid
  julia --project=perf perf/flame.jl --method cube --scenario hex_r5_d50
  julia --project=perf perf/flame.jl --resolution 8 --n-iterations 20
=#

include("common.jl")
import Profile
import ProfileCanvas

function parse_flame_args()
    args = parse_perf_args(; description="Generate flame graph for gridPathCells")
    return args
end

function profile_method!(scenarios, n_iterations, method_fn, label)
    for sc in scenarios
        method_fn(sc.start_cell, sc.end_cell)
    end

    Profile.init(; delay=0.00001)
    Profile.clear()
    Profile.@profile for _ in 1:n_iterations
        for sc in scenarios
            method_fn(sc.start_cell, sc.end_cell)
        end
    end

    data = Profile.fetch()
    if isempty(data)
        println("Warning: No profile samples collected for $label. Try increasing --n-iterations.")
        return
    end

    output_dir = joinpath(@__DIR__, "output")
    mkpath(output_dir)
    outfile = joinpath(output_dir, "flame_$(label).html")
    ProfileCanvas.html_file(outfile, data)
    println("Flame graph written to: $outfile")
end

function main()
    args = parse_flame_args()
    n_iterations = args["n-iterations"]
    scenarios = filter_scenarios(build_scenarios(), args)

    if isempty(scenarios)
        println("No scenarios matched the filters.")
        return
    end

    method = get(args, "method", "all")
    Ext = FastH3.FastH3Extension

    suffix = args["scenario"]
    tag = isempty(suffix) ? "" : "_$(suffix)"

    if method in ("cube", "all")
        println("Profiling cube gridPathCells ($(length(scenarios)) scenarios, $(n_iterations) iterations)")
        cube_fn = (s, e) -> FastH3.gridPathCells(s, e)
        profile_method!(scenarios, n_iterations, cube_fn, "gridPathCells_cube$(tag)")
    end

    if method in ("gc", "all")
        println("Profiling great-circle gridPathCells ($(length(scenarios)) scenarios, $(n_iterations) iterations)")
        gc_fn = (s, e) -> Ext.gridPathCells(s, e)
        profile_method!(scenarios, n_iterations, gc_fn, "gridPathCells_gc$(tag)")
    end

    if method in ("hybrid", "all")
        println("Profiling hybrid gridPathCellsHybrid ($(length(scenarios)) scenarios, $(n_iterations) iterations)")
        hyb_fn = (s, e) -> Ext.gridPathCellsHybrid(s, e)
        profile_method!(scenarios, n_iterations, hyb_fn, "gridPathCells_hybrid$(tag)")
    end
end

main()
