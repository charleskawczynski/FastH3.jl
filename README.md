# FastH3.jl

| **Build Status**                                                                                                                                          | **Documentation**                                                                                                                                                                                                                          | **Coverage**                                                                                                                                  |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------:|
| [![CI](https://github.com/charleskawczynski/FastH3.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/charleskawczynski/FastH3.jl/actions/workflows/ci.yml) | [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://charleskawczynski.github.io/FastH3.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://charleskawczynski.github.io/FastH3.jl/dev/) | [![Codecov](https://codecov.io/gh/charleskawczynski/FastH3.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/charleskawczynski/FastH3.jl) |

A pure Julia implementation of [Uber's H3](https://github.com/uber/h3) hexagonal hierarchical geospatial indexing system.

## Acknowledgments

This library is a direct translation of the [H3 C library](https://github.com/uber/h3) developed by Uber. The algorithms, lookup tables, and index layout are faithfully ported from the C reference implementation.

[H3.jl](https://github.com/wookay/H3.jl) is the existing Julia package for H3. It wraps the H3 C library via `ccall` and `BinaryProvider`, providing Julia access to the original C implementation.

## Why FastH3.jl?

FastH3.jl is a pure Julia reimplementation with zero binary or C dependencies. This matters for several reasons:

- **No binary dependencies** -- H3.jl depends on `BinaryProvider` to distribute pre-compiled `libh3` binaries. FastH3.jl has no compiled C code, no JLL packages, and no platform-specific binaries. This eliminates build failures on unsupported or exotic platforms and simplifies installation.
- **Julia-native composability** -- Pure Julia code is transparent to the compiler. Functions can be specialized via multiple dispatch, inlined across call boundaries, and composed with Julia's type system (e.g., custom number types, autodiff) in ways that are impossible through a `ccall` FFI boundary.
- **GPU / accelerator portability** -- Pure Julia code can, in principle, be compiled to GPU kernels (via CUDA.jl, AMDGPU.jl, KernelAbstractions.jl) or run on other accelerators. C-bound code cannot cross that boundary.
- **No heap-allocating `Ref` wrappers** -- The H3 C API writes results through output pointers, so H3.jl must allocate a mutable `Ref{T}()` (a heap-allocated `Base.RefValue{T}`) for every call that returns a scalar -- `Ref{H3Index}()`, `Ref{Int64}()`, `Ref{Cdouble}()`, `Ref{LatLng}()`, `Ref{CellBoundary}()`, etc. FastH3.jl returns results directly as tuples, avoiding those per-call heap allocations entirely.
- **Auditability and hackability** -- All logic lives in Julia source files, making it straightforward to read, modify, debug with Julia's tooling, and contribute to without a C toolchain.
