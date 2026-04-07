# Callback enumeration (no result vectors)

Several APIs can return many cells (or related indexes) in a `Vector`. For each of these,
FastH3 also provides an overload whose **first** argument is a callback `op`. That overload
returns only `H3Error` and invokes `op` once per output element in the same order as the
vector-returning API would produce. The callback path does not allocate the result
collection itself (your `op` may still allocate).

Use Julia’s `do` syntax so the callback reads as a trailing block:

```julia
using FastH3

err, h = latLngToCell(LatLng(0.0, 0.0), 8)
err == E_SUCCESS || error("latLngToCell failed")

# Collect children without a dedicated cellToChildren vector allocation
# (still allocates the collector `out` if you push!)
out = H3Index[]
err = cellToChildren(h, 9) do child
    push!(out, child)
end
err == E_SUCCESS || error("cellToChildren failed")
```

## APIs with callback overloads

| Callback form | Typical use |
|---------------|-------------|
| `cellToChildren(op, h, childRes)` | All children at finer resolution |
| `uncompactCells(op, compactedSet, res)` | Expand compacted set to `res` |
| `getPentagons(op, res)` | Pentagon base cells at `res` |
| `getRes0Cells(op)` | All resolution-0 indexes |
| `getIcosahedronFaces(op, h)` | Face indices (see note below) |
| `gridDiskDistancesUnsafe(op, origin, k)` | `op(cell, ring)` for disk + distance |
| `gridDiskUnsafe(op, origin, k)` | `op(cell)` only |
| `gridRingUnsafe(op, origin, k)` | Ring at exact distance `k` |
| `gridPathCells(op, start, end)` | Cells along interpolated path |
| `cellToVertexes(op, origin)` | Vertex indexes for a cell |
| `originToDirectedEdges(op, origin)` | Directed edges from a cell |

Docstrings for these functions describe the callback signature and link to the H3 C API.

## `getIcosahedronFaces` and the vector form

The `Vector{Int}` overload returns a fixed-length vector aligned with `maxFaceCount(h)`, using
`-1` for unused slots (internal sentinel). The **callback** overload calls `op(face)` only for
**valid** face indices (sentinel slots are skipped). If you need the raw buffer including
sentinels, use the vector-returning method.

## What does not have a callback overload yet

`gridDisk` and `gridDiskDistances` may fall back to an internal BFS when the fast “unsafe” path
hits a pentagon. That path uses extra heap structures, so a zero-allocation callback story is
not offered for those entry points yet. For allocation-free disk enumeration when the unsafe
algorithm applies, use `gridDiskUnsafe` or `gridDiskDistancesUnsafe` with a callback.

`compactCells`, polygon-to-cells routines, and similar algorithms are also callback-free for now.

## See also

- [Hierarchy](hierarchy.md) — parent/child and compact/uncompact
- [Traversal](traversal.md) — grid disk, ring, paths
- [Inspection](inspection.md) — `getIcosahedronFaces`, `getRes0Cells`, and related helpers
