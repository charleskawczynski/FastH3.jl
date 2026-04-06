# Extensions

`FastH3.FastH3Extension` provides alternative `gridPathCells` implementations
that complement the core H3 cube-interpolation algorithm. All functions live in
the `FastH3Extension` submodule and must be accessed with a qualified name:

```julia
import FastH3
path = FastH3.FastH3Extension.gridPathCells(c0, c1)
```

## Algorithm Comparison

FastH3 ships three ways to compute a path of H3 cells between two endpoints.
They differ in algorithm, robustness, speed, and the semantic meaning of the
path they return.

### `FastH3.gridPathCells` — Cube Interpolation

The standard H3 algorithm. It maps both endpoints into a shared local IJ
coordinate system, computes `gridDistance`, and linearly interpolates in cube
coordinates. Each interpolated point is rounded to the nearest cell.

The path follows a straight line in IJK coordinate space, which does not
correspond to the geographic great circle on the sphere:

```mermaid
graph LR
    S["Start ●"] --> A["·"] --> B["·"] --> C["·"] --> E["End ●"]
```

The cube path takes the shortest route through the hex grid (minimum hops).
Cells that the great-circle arc would cross but that lie off the IJK line
are not included.

| Property | Value |
|:--|:--|
| **Path length** | Exactly `gridDistance + 1` cells (the shortest grid path). |
| **Path meaning** | Follows the IJK coordinate grid, not the geographic great circle. |
| **Speed** | Very fast — a single allocation and O(n) integer arithmetic. |
| **Robustness** | Returns `E_DOMAIN` when start and end lie on different icosahedral faces that cannot share a local IJ coordinate system. Can also silently produce incorrect endpoints for certain cross-face pairs where `gridDistance` succeeds but the IJ mapping is inaccurate. |

```julia
err, path = FastH3.gridPathCells(c0, c1)
# err may be E_DOMAIN for cross-face paths
```

### `FastH3.FastH3Extension.gridPathCells` — Great-Circle Walk

An analytic walk along the geographic great-circle arc connecting the two cell
centres. At each cell it computes the intersection of the great-circle arc with
every boundary edge, steps into the neighbor on the far side of the nearest
exit edge, and repeats. A bisection fallback handles degenerate cases (e.g. the
arc grazing a hex vertex).

The path follows the great-circle arc on the sphere and includes every cell
whose boundary it crosses, even corner cells that the shortest grid path would
skip:

```mermaid
graph LR
    S["Start ●"] --> A["·"] --> X["· (corner clip)"] --> B["·"] --> C["·"] --> E["End ●"]
```

The step-by-step walk works as follows:

```mermaid
flowchart TD
    fetchBoundary["Fetch current cell boundary vertices"] --> intersect["Compute arc ∩ each boundary edge"]
    intersect --> pickExit["Pick nearest exit crossing past current position"]
    pickExit --> nudge["Nudge exit point outward → latLngToCell"]
    nudge --> nextCell["Step into neighbor cell"]
    nextCell --> check{"Reached destination?"}
    check -->|No| fetchBoundary
    check -->|Yes| done["Done"]
```

| Property | Value |
|:--|:--|
| **Path length** | Variable — includes every cell whose boundary the arc crosses. May be longer than `gridDistance + 1` because the arc can clip corner cells that the shortest grid path would skip. |
| **Path meaning** | Geographically faithful: the path follows the great-circle arc on the sphere. |
| **Speed** | Slower — each step requires trigonometric intersection tests against up to 6 boundary edges plus `latLngToCell` lookups. Roughly 10–25x slower than cube interpolation. |
| **Robustness** | Always succeeds, including across icosahedral face boundaries and near pentagons. |

```julia
path = FastH3.FastH3Extension.gridPathCells(c0, c1)
# always returns a valid path
```

A callback form is also available for zero-allocation iteration:

```julia
FastH3.FastH3Extension.gridPathCells!(c0, c1) do cell
    # process each H3Index
end
```

### `FastH3.FastH3Extension.gridPathCellsHybrid` — Piecewise Cube with GC Fallback

A hybrid that gets cube speed for the common case and GC robustness for
cross-face paths.

The decision logic at each recursion level:

```mermaid
flowchart TD
    entry["gridPathCellsHybrid(c0, c1)"] --> tryCube{"cube(c0, c1) succeeds\nand ends at c1?"}
    tryCube -->|Yes| returnCube["Return cube path"]
    tryCube -->|No| trivial{"c0 == c1\nor neighbors?"}
    trivial -->|Yes| returnDirect["Return directly"]
    trivial -->|No| slerp["mid = SLERP midpoint → latLngToCell"]
    slerp --> checkMid{"mid == c0\nor mid == c1?"}
    checkMid -->|Yes| fallbackGC["Fall back to GC walk"]
    checkMid -->|No| recurseLeft["hybrid(c0, mid)"]
    recurseLeft --> recurseRight["hybrid(mid, c1)"]
    recurseRight --> stitch["Stitch: left ++ right"]
```

For cross-face paths, the recursion splits at the great-circle midpoint and
applies cube interpolation to each same-face sub-segment:

```mermaid
graph LR
    subgraph faceA ["Face A — cube interpolation"]
        S["Start ●"] --> a1["·"] --> a2["·"] --> a3["·"]
    end
    a3 --> mid["Mid (SLERP)"]
    subgraph faceB ["Face B — cube interpolation"]
        mid --> b1["·"] --> b2["·"] --> E["End ●"]
    end
```

The recursion depth is O(log n) where n is the number of cells, bounded by
`max_depth`. Each face-local segment uses fast cube interpolation. Only the
handful of cells at face boundaries may trigger a GC walk fallback.

| Property | Value |
|:--|:--|
| **Path length** | Piecewise-shortest within each icosahedral face. The total may differ from both the global `gridDistance` and the great-circle cell count. |
| **Path meaning** | Within each face, follows cube (IJK) interpolation. Face-boundary waypoints are determined by the great circle. Neither the global shortest grid path nor the great-circle-faithful path. |
| **Speed** | Same-face: identical to cube. Cross-face: ~2–8x faster than the pure great-circle walk (O(log n) bisections, then cube for each sub-segment). |
| **Robustness** | Always succeeds — inherits GC walk robustness for any segment that cube cannot handle. |

```julia
path = FastH3.FastH3Extension.gridPathCellsHybrid(c0, c1)
# always returns a valid path
```

## Summary Table

| | Cube | Great-Circle Walk | Hybrid |
|:--|:--|:--|:--|
| **Function** | `FastH3.gridPathCells` | `FastH3.FastH3Extension.gridPathCells` | `FastH3.FastH3Extension.gridPathCellsHybrid` |
| **Returns** | `(H3Error, Vector{H3Index})` | `Vector{H3Index}` | `Vector{H3Index}` |
| **Always succeeds** | No | Yes | Yes |
| **Shortest path** | Yes (when it works) | No | Piecewise per face |
| **Follows great circle** | No | Yes | No |
| **Cross-face** | Fails or wrong result | Works | Works |
| **Relative speed** | 1x | 10–25x slower | ~1x same-face, ~2–8x slower cross-face |

## When to Use Which

- **Same-face, performance-critical**: Use `FastH3.gridPathCells`. It is the
  fastest option and produces the shortest grid path, but only when both cells
  share an icosahedral face.
- **Need the geographic great-circle path**: Use
  `FastH3.FastH3Extension.gridPathCells`. It faithfully traces the spherical
  arc, which matters for applications like line-of-sight, coverage analysis, or
  trajectory mapping where the physical path on the globe is significant.
- **General-purpose robust pathing**: Use
  `FastH3.FastH3Extension.gridPathCellsHybrid`. It automatically uses the fast
  cube algorithm where possible and gracefully handles cross-face paths. A good
  default when you need reliability without paying full great-circle cost.

## API Reference

```@docs
FastH3.FastH3Extension.gridPathCells!
FastH3.FastH3Extension.gridPathCells
FastH3.FastH3Extension.gridPathCellsHybrid
```
