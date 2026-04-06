#=
Shared utilities for FastH3 performance scripts.

Generates a wide range of (start, end) cell pairs for benchmarking
gridPathCells across resolutions, distances, pentagon cells, and
cross-base-cell paths.

Common CLI flags (parsed by `parse_perf_args`):
  --resolution   H3 resolution to benchmark (default: all)
  --distance     Grid distance to benchmark (default: all)
  --n-iterations Iteration count for profiling (default: 5)
  --scenario     Scenario filter substring (default: run all)
=#

import ArgParse
import FastH3

struct Scenario
    name::String
    start_cell::FastH3.H3Index
    end_cell::FastH3.H3Index
    distance::Int
    resolution::Int
end

function _cell_at_distance(origin::FastH3.H3Index, k::Int)::Union{FastH3.H3Index, Nothing}
    k == 0 && return origin
    err, ring = FastH3.gridRingUnsafe(origin, k)
    if err == FastH3.E_SUCCESS && !isempty(ring)
        for cell in ring
            cell != FastH3.H3_NULL && return cell
        end
    end
    err_d, disk = FastH3.gridDisk(origin, k)
    err_d != FastH3.E_SUCCESS && return nothing
    for cell in disk
        cell == FastH3.H3_NULL && continue
        err_g, dist = FastH3.gridDistance(origin, cell)
        err_g == FastH3.E_SUCCESS && dist == k && return cell
    end
    return nothing
end

const BENCH_ORIGIN = FastH3.LatLng(deg2rad(40.7128), deg2rad(-74.0060))  # NYC
const BENCH_RESOLUTIONS = [0, 3, 5, 8, 10, 15]
const BENCH_DISTANCES = [0, 1, 2, 5, 10, 20, 50]

function build_scenarios(;
    resolutions=BENCH_RESOLUTIONS,
    distances=BENCH_DISTANCES,
    include_pentagon::Bool=true,
    include_cross_base_cell::Bool=true,
)::Vector{Scenario}
    scenarios = Scenario[]

    for res in resolutions
        err, origin = FastH3.latLngToCell(BENCH_ORIGIN, res)
        err != FastH3.E_SUCCESS && continue

        for dist in distances
            target = _cell_at_distance(origin, dist)
            target === nothing && continue
            push!(scenarios, Scenario("hex_r$(res)_d$(dist)", origin, target, dist, res))
        end
    end

    if include_pentagon
        for res in resolutions
            err, pentagons = FastH3.getPentagons(res)
            err != FastH3.E_SUCCESS && continue
            pent = pentagons[1]
            for dist in [0, 1, 2, 3]
                target = _cell_at_distance(pent, dist)
                target === nothing && continue
                push!(scenarios, Scenario("pent_r$(res)_d$(dist)", pent, target, dist, res))
            end
        end
    end

    if include_cross_base_cell
        for res in [3, 5, 8]
            err, origin = FastH3.latLngToCell(BENCH_ORIGIN, res)
            err != FastH3.E_SUCCESS && continue
            originBC = FastH3.h3_get_base_cell(origin)

            max_k = res <= 3 ? 15 : 10
            err_d, disk = FastH3.gridDisk(origin, max_k)
            err_d != FastH3.E_SUCCESS && continue
            found = Set{Int}()
            for cell in disk
                cell == FastH3.H3_NULL && continue
                FastH3.h3_get_base_cell(cell) == originBC && continue
                err_g, dist = FastH3.gridDistance(origin, cell)
                err_g != FastH3.E_SUCCESS && continue
                dist in found && continue
                push!(found, dist)
                push!(scenarios, Scenario("xbase_r$(res)_d$(dist)", origin, cell, dist, res))
            end
        end
    end

    return scenarios
end

function parse_perf_args(; description="FastH3 performance script")
    s = ArgParse.ArgParseSettings(; description)
    ArgParse.@add_arg_table! s begin
        "--resolution"
        help = "Single H3 resolution to benchmark (omit for all)"
        arg_type = Int
        default = -1
        "--distance"
        help = "Single grid distance to benchmark (omit for all)"
        arg_type = Int
        default = -1
        "--n-iterations"
        help = "Number of iterations for profiling"
        arg_type = Int
        default = 10_000
        "--scenario"
        help = "Substring filter on scenario name"
        default = ""
        "--method"
        help = "Which implementation to profile: cube, gc, robust, or all"
        default = "all"
    end
    return ArgParse.parse_args(s)
end

function filter_scenarios(scenarios::Vector{Scenario}, args)
    out = scenarios
    res = args["resolution"]
    dist = args["distance"]
    substr = args["scenario"]

    if res >= 0
        out = filter(s -> s.resolution == res, out)
    end
    if dist >= 0
        out = filter(s -> s.distance == dist, out)
    end
    if !isempty(substr)
        out = filter(s -> occursin(substr, s.name), out)
    end
    return out
end
