# H3X.jl

*A pure Julia implementation of Uber's [H3](https://h3geo.org/) hexagonal hierarchical geospatial indexing system.*

## Overview

H3X.jl provides a complete, dependency-free Julia implementation of the H3 library.
It maps latitude/longitude coordinates onto a hierarchical grid of hexagonal (and
pentagonal) cells, enabling efficient spatial indexing, neighbor traversal, and
distance computation.

Key properties:

- **Pure Julia** — no C dependencies or binary wrappers.
- **Zero-allocation scalars** — core index operations allocate nothing on the heap.
- **Full H3 v4 API** — covers indexing, inspection, traversal, hierarchy, directed
  edges, vertexes, and measurement functions.

## Quick Start

```julia
using H3X

lat, lng = 0.6518070561696664, -1.7453292519943295  # radians
err, cell = latLngToCell(LatLng(lat, lng), 5)
err, center = cellToLatLng(cell)
err, disk = gridDisk(cell, 1)
```

## C API Correspondence

Every public function in H3X.jl corresponds to a function in the
[H3 C library](https://h3geo.org/docs/api/indexing). Docstrings include direct
links to the matching C API documentation page.

## Module

```@docs
H3X.H3X
```

## Contents

```@contents
Pages = [
    "api/indexing.md",
    "api/inspection.md",
    "api/traversal.md",
    "api/hierarchy.md",
    "api/directed_edges.md",
    "api/vertexes.md",
    "api/measurement.md",
]
Depth = 2
```
