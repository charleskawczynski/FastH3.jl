#=
Great-Circle Path — analytic great-circle walk for gridPathCells
================================================================

Given two H3 cells, compute the sequence of H3 cells that a great-circle
arc between their centres traverses.

Algorithm
---------
The path is a great-circle arc on the unit sphere.  Each H3 cell boundary
is a polygon whose edges are also great-circle arcs.  Two great circles
intersect at exactly two antipodal points, obtainable in closed form from
the cross product of their normal vectors:

    n_path = u0 × u1          (normal to the path's great circle)
    n_edge = vᵢ × vⱼ          (normal to a hex-edge great circle)
    I      = ± normalize(n_path × n_edge)

We walk the path cell-by-cell:

  1. Fetch the current cell's boundary vertices (in Cartesian).
  2. For each edge, compute the two candidate intersection points.
  3. Keep only candidates that lie on *both* the path arc and the edge
     arc (checked via the `_on_arc` predicate).
  4. Among valid crossings past the current position, pick the one with
     the smallest arc-length fraction — that is the exit point.
  5. Step into the neighbor on the far side of the crossed edge.
  6. Repeat until we reach the destination cell.

Neighbor resolution (step 5) nudges the exit point outward from the
current cell's centroid by 1 % of the edge's angular length and calls
`latLngToCell`.

If the analytic step fails (e.g. the path grazes a hex vertex where
three cells meet), a bisection fallback handles the remainder.
=#

# ── Vec3 type and inline helpers ──

const Vec3 = NTuple{3,Float64}

@inline _cross3(a::Vec3, b::Vec3)::Vec3 =
    (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])

@inline _dot3(a::Vec3, b::Vec3)::Float64 = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _norm3(a::Vec3)::Float64 = sqrt(_dot3(a, a))
@inline _scale3(s::Float64, a::Vec3)::Vec3 = (s * a[1], s * a[2], s * a[3])
@inline _neg3(a::Vec3)::Vec3 = (-a[1], -a[2], -a[3])

# ── Geo ↔ Cartesian conversions ──

function _latLngToVec3(ll::LatLng)::Vec3
    clat = cos(ll.lat)
    return (clat * cos(ll.lng), clat * sin(ll.lng), sin(ll.lat))
end

function _vec3ToLatLng(v::Vec3)::LatLng
    x, y, z = v
    return LatLng(atan(z, sqrt(x * x + y * y)), atan(y, x))
end

function _cellToBoundaryVec3(cell::H3Index)::Tuple{Int,NTuple{10,Vec3}}
    _, cb = cellToBoundary(cell)
    nv = Int(cb.numVerts)
    cart_verts = ntuple(Val(10)) do i
        _latLngToVec3(cb.verts[i])
    end
    return (nv, cart_verts)
end

# ── SLERP interpolation along the great-circle arc ──

function _interpolate_great_circle(p0::LatLng, p1::LatLng, frac::Float64)::LatLng
    frac <= 0.0 && return p0
    frac >= 1.0 && return p1
    u0 = _latLngToVec3(p0)
    u1 = _latLngToVec3(p1)
    cos_omega = clamp(_dot3(u0, u1), -1.0, 1.0)
    omega = acos(cos_omega)
    omega < 1e-15 && return p0
    sin_omega = sin(omega)
    a = sin((1.0 - frac) * omega) / sin_omega
    b = sin(frac * omega) / sin_omega
    v = (a * u0[1] + b * u1[1], a * u0[2] + b * u1[2], a * u0[3] + b * u1[3])
    return _vec3ToLatLng(v)
end

# ── Analytic great-circle intersection helpers ──

const _ARC_EPS = 1e-12
const _FRAC_EPS = 1e-10

"""
P lies on the shorter great-circle arc from A to B.
`n_AB = _cross3(A, B)` must be precomputed.
"""
@inline function _on_arc(P::Vec3, A::Vec3, B::Vec3, n_AB::Vec3)::Bool
    return _dot3(_cross3(A, P), n_AB) >= -_ARC_EPS &&
           _dot3(_cross3(B, P), n_AB) <= _ARC_EPS
end

"""
Fraction ∈ [0,1] of the great-circle arc u0→u1 where P lies.
"""
@inline function _arc_fraction(u0::Vec3, u1::Vec3, P::Vec3)::Float64
    d_total = acos(clamp(_dot3(u0, u1), -1.0, 1.0))
    d_total < 1e-15 && return 0.0
    d_0P = acos(clamp(_dot3(u0, P), -1.0, 1.0))
    return d_0P / d_total
end

"""
Find where the great-circle path u0→u1 exits through a cell boundary.

Returns `(found, exit_frac, exit_pt, edge_vi, edge_vj)`.
"""
function _find_exit_crossing(
    u0::Vec3, u1::Vec3, n_path::Vec3,
    nv::Int, verts::NTuple{10,Vec3}, min_frac::Float64,
)::Tuple{Bool,Float64,Vec3,Vec3,Vec3}
    best_frac = 2.0
    best_pt = u1
    best_vi = u0
    best_vj = u1
    found = false

    for i in 1:nv
        vi = verts[i]
        vj = verts[mod1(i + 1, nv)]
        n_edge = _cross3(vi, vj)

        ix = _cross3(n_path, n_edge)
        n_ix = _norm3(ix)
        n_ix < _ARC_EPS && continue

        inv_n = 1.0 / n_ix
        candidate = _scale3(inv_n, ix)

        for sign in (1.0, -1.0)
            pt = sign > 0 ? candidate : _neg3(candidate)

            _on_arc(pt, vi, vj, n_edge) || continue
            _on_arc(pt, u0, u1, n_path) || continue

            frac = _arc_fraction(u0, u1, pt)
            frac <= min_frac + _FRAC_EPS && continue

            if frac < best_frac
                best_frac = frac
                best_pt = pt
                best_vi = vi
                best_vj = vj
                found = true
            end
        end
    end

    return (found, best_frac, best_pt, best_vi, best_vj)
end

"""
Determine the H3 cell on the far side of the crossed edge.

Nudges `exit_pt` outward from the current cell's centroid by 1 % of the
edge's angular length, then calls `latLngToCell`.
"""
function _next_cell_across_edge(
    exit_pt::Vec3, vi::Vec3, vj::Vec3,
    nv::Int, cell_verts::NTuple{10,Vec3}, h3_res::Int,
)::H3Index
    cx, cy, cz = 0.0, 0.0, 0.0
    for i in 1:nv
        v = cell_verts[i]
        cx += v[1]
        cy += v[2]
        cz += v[3]
    end
    centroid = (cx / nv, cy / nv, cz / nv)

    dx = exit_pt[1] - centroid[1]
    dy = exit_pt[2] - centroid[2]
    dz = exit_pt[3] - centroid[3]
    d_norm = sqrt(dx * dx + dy * dy + dz * dz)

    if d_norm < 1e-15
        _, cell = latLngToCell(_vec3ToLatLng(exit_pt), h3_res)
        return cell
    end

    edge_angle = acos(clamp(_dot3(vi, vj), -1.0, 1.0))
    nudge_scale = edge_angle * 0.01 / d_norm

    nx = exit_pt[1] + nudge_scale * dx
    ny = exit_pt[2] + nudge_scale * dy
    nz = exit_pt[3] + nudge_scale * dz
    inv_n = 1.0 / sqrt(nx * nx + ny * ny + nz * nz)

    _, cell = latLngToCell(_vec3ToLatLng((nx * inv_n, ny * inv_n, nz * inv_n)), h3_res)
    return cell
end

# ── Bisection fallback ──

function _no_cells_between(
    h3_res::Int,
    frac0::Float64, c0::H3Index,
    frac1::Float64, c1::H3Index,
    p_start::LatLng, p_end::LatLng,
    n_iter::Int, max_iter::Int,
)::Bool
    n_iter >= max_iter && return true
    frac_mid = (frac0 + frac1) * 0.5
    p_mid = _interpolate_great_circle(p_start, p_end, frac_mid)
    _, cm = latLngToCell(p_mid, h3_res)
    (cm != c0 && cm != c1) && return false
    _no_cells_between(h3_res, frac0, c0, frac_mid, cm, p_start, p_end, n_iter + 1, max_iter) || return false
    return _no_cells_between(h3_res, frac_mid, cm, frac1, c1, p_start, p_end, n_iter + 1, max_iter)
end

function _find_edge_crossing_frac(
    h3_res::Int,
    frac0::Float64, c0::H3Index,
    frac1::Float64, c1::H3Index,
    p_start::LatLng, p_end::LatLng,
    n_iter::Int, max_iter::Int,
)::Float64
    frac_mid = (frac0 + frac1) * 0.5
    n_iter >= max_iter && return frac_mid
    p_mid = _interpolate_great_circle(p_start, p_end, frac_mid)
    _, cm = latLngToCell(p_mid, h3_res)
    cm == c0 && return _find_edge_crossing_frac(h3_res, frac_mid, cm, frac1, c1, p_start, p_end, n_iter + 1, max_iter)
    cm == c1 && return _find_edge_crossing_frac(h3_res, frac0, c0, frac_mid, cm, p_start, p_end, n_iter + 1, max_iter)
    return frac_mid
end

function _grid_path_fallback!(
    emit!,
    h3_res::Int,
    frac0::Float64, c0::H3Index,
    frac1::Float64, c1::H3Index,
    p_start::LatLng, p_end::LatLng,
    depth::Int, max_depth::Int, max_iter::Int,
)
    if c0 == c1
        emit!(c0)
        return nothing
    end

    _, are_neighbors = areNeighborCells(c0, c1)
    if are_neighbors && _no_cells_between(h3_res, frac0, c0, frac1, c1, p_start, p_end, 0, max_iter)
        emit!(c0)
        emit!(c1)
        return nothing
    end

    if depth <= max_depth
        frac_mid = (frac0 + frac1) * 0.5
        p_mid = _interpolate_great_circle(p_start, p_end, frac_mid)
        _, cm = latLngToCell(p_mid, h3_res)
        _grid_path_fallback!(emit!, h3_res, frac0, c0, frac_mid, cm, p_start, p_end, depth + 1, max_depth, max_iter)
        _grid_path_fallback!(emit!, h3_res, frac_mid, cm, frac1, c1, p_start, p_end, depth + 1, max_depth, max_iter)
    end
    return nothing
end

# ── Neighbor validation ──

"""
Check if `c1` is a geographic neighbor of `c0`.

Uses `areNeighborCells` first, then falls back to a geographic proximity
test for cross-icosahedral-face boundaries where H3's local-IJ coordinate
system cannot connect the cells.  Two cells whose centres are within
`2 √3` average hex edge lengths are treated as neighbors.
"""
function _is_geographic_neighbor(c0::H3Index, c1::H3Index)::Bool
    _, is_n = areNeighborCells(c0, c1)
    is_n && return true
    _, p0 = cellToLatLng(c0)
    _, p1 = cellToLatLng(c1)
    dist_rads = greatCircleDistanceRads(p0, p1)
    dist_km = dist_rads * EARTH_RADIUS_KM
    res = getResolution(c0)
    _, avg_edge_km = getHexagonEdgeLengthAvgKm(res)
    return dist_km < avg_edge_km * 2.0
end

# ── Arc-based nudge for finding the next cell ──

"""
When the centroid-outward nudge resolves to a non-neighbor (common at
icosahedral face boundaries), sample points slightly ahead on the great-circle
arc to find the correct next cell.
"""
function _try_arc_nudge(
    u0::Vec3, u1::Vec3, exit_frac::Float64,
    current_cell::H3Index, h3_res::Int,
)::H3Index
    omega = acos(clamp(_dot3(u0, u1), -1.0, 1.0))
    omega < 1e-15 && return current_cell
    sin_omega = sin(omega)

    for nudge in (1e-8, 1e-6, 1e-4, 1e-3, 5e-3)
        t = min(exit_frac + nudge, 1.0)
        a = sin((1.0 - t) * omega) / sin_omega
        b = sin(t * omega) / sin_omega
        pt = (a * u0[1] + b * u1[1], a * u0[2] + b * u1[2], a * u0[3] + b * u1[3])
        _, cell = latLngToCell(_vec3ToLatLng(pt), h3_res)
        cell == current_cell && continue
        _is_geographic_neighbor(current_cell, cell) && return cell
    end
    return current_cell
end

# ── Main great-circle walk ──

function _grid_path_walk!(
    emit!, start_::H3Index, end_::H3Index,
    max_depth::Int, max_iter::Int,
)
    h3_res = getResolution(start_)
    _, p0 = cellToLatLng(start_)
    _, p1 = cellToLatLng(end_)
    c1 = end_

    if start_ == c1
        emit!(start_)
        return nothing
    end

    u0 = _latLngToVec3(p0)
    u1 = _latLngToVec3(p1)
    n_path = _cross3(u0, u1)

    if _norm3(n_path) < 1e-15
        _grid_path_fallback!(emit!, h3_res, 0.0, start_, 1.0, c1, p0, p1, 0, max_depth, max_iter)
        return nothing
    end

    current_cell = start_
    current_frac = 0.0

    for _ in 1:max_depth
        current_cell == c1 && break

        (nv, verts) = _cellToBoundaryVec3(current_cell)
        (found, exit_frac, exit_pt, edge_vi, edge_vj) =
            _find_exit_crossing(u0, u1, n_path, nv, verts, current_frac)

        if !found || exit_frac > 1.0 + _FRAC_EPS
            _grid_path_fallback!(emit!, h3_res, current_frac, current_cell, 1.0, c1, p0, p1, 0, max_depth, max_iter)
            return nothing
        end

        exit_frac = min(exit_frac, 1.0)
        emit!(current_cell)

        next_cell = _next_cell_across_edge(exit_pt, edge_vi, edge_vj, nv, verts, h3_res)

        if next_cell == current_cell
            arc_cell = _try_arc_nudge(u0, u1, exit_frac, current_cell, h3_res)
            if arc_cell != current_cell
                next_cell = arc_cell
            else
                _grid_path_fallback!(emit!, h3_res, exit_frac, current_cell, 1.0, c1, p0, p1, 0, max_depth, max_iter)
                return nothing
            end
        end

        current_cell = next_cell
        current_frac = exit_frac
    end

    emit!(current_cell)
    return nothing
end

# ── Gap-filling post-processing ──

"""
Fill a non-neighbor gap between `c0` and `c1` by bisecting the great-circle
arc between their centres and recursing.  Emits intermediate cells via `f`.
"""
function _fill_gap!(f, c0::H3Index, c1::H3Index, h3_res::Int, depth::Int)
    depth > 50 && return
    _is_geographic_neighbor(c0, c1) && return
    c0 == c1 && return

    _, p0 = cellToLatLng(c0)
    _, p1 = cellToLatLng(c1)

    for frac in (0.5, 0.33, 0.67, 0.25, 0.75)
        pm = _interpolate_great_circle(p0, p1, frac)
        _, cm = latLngToCell(pm, h3_res)
        if cm != c0 && cm != c1
            _fill_gap!(f, c0, cm, h3_res, depth + 1)
            f(cm)
            _fill_gap!(f, cm, c1, h3_res, depth + 1)
            return
        end
    end
end

# ── Public API: Great-Circle Walk ──

"""
    gridPathCells!(f, start_::H3Index, end_::H3Index; max_depth=100, max_iter=10)

Walk the great-circle arc between the centres of `start_` and `end_` and call
`f(cell::H3Index)` for each H3 cell traversed, in order.  The first call
receives `start_` and the last receives `end_`.

Uses an analytic great-circle / cell-boundary intersection walk with a
bisection fallback for degenerate cases.

`max_depth` limits the walk and fallback recursion depth.
`max_iter` limits iterations inside the bisection neighbour check.
"""
function gridPathCells!(f, start_::H3Index, end_::H3Index; max_depth::Int=100, max_iter::Int=10)
    h3_res = getResolution(start_)
    prev = Ref(H3_NULL)
    _grid_path_walk!(start_, end_, max_depth, max_iter) do cell
        if cell != prev[]
            if prev[] != H3_NULL && !_is_geographic_neighbor(prev[], cell)
                _fill_gap!(f, prev[], cell, h3_res, 0)
            end
            f(cell)
            prev[] = cell
        end
    end
    return nothing
end

"""
    gridPathCells(start_::H3Index, end_::H3Index; max_depth=100, max_iter=10) -> Vector{H3Index}

Return the ordered sequence of H3 cells that the great-circle arc between the
centres of `start_` and `end_` passes through.  Equivalent to collecting the
cells emitted by [`gridPathCells!`](@ref).
"""
function gridPathCells(start_::H3Index, end_::H3Index; max_depth::Int=100, max_iter::Int=10)
    path = H3Index[]
    gridPathCells!(start_, end_; max_depth, max_iter) do cell
        push!(path, cell)
    end
    return path
end

# ── Public API: Robust (cube with greedy re-anchoring) ──

"""
Pick the neighbor of `current` that is geographically closest to `target`.
Uses `gridDisk(current, 1)` to enumerate the 6 (or 5) immediate neighbors,
then selects the one minimising great-circle distance to `target`.
"""
function _step_toward(current::H3Index, target::H3Index)::H3Index
    _, target_ll = cellToLatLng(target)
    err, disk = gridDisk(current, 1)
    err != E_SUCCESS && return current
    best = current
    best_dist = Inf
    for cell in disk
        cell == H3_NULL && continue
        cell == current && continue
        _, cell_ll = cellToLatLng(cell)
        d = greatCircleDistanceRads(cell_ll, target_ll)
        if d < best_dist
            best_dist = d
            best = cell
        end
    end
    return best
end

"""
    gridPathCellsRobust(start_::H3Index, end_::H3Index) -> Vector{H3Index}

Compute a path of H3 cells between `start_` and `end_` using cube
interpolation with greedy re-anchoring at icosahedral face boundaries.

Tries `FastH3.gridPathCells` (cube) from the current position toward `end_`.
When cube fails (cross-face `E_DOMAIN` or incorrect endpoint due to IJ mapping
errors), advances one cell toward `end_` by picking the geographically closest
neighbor, then retries cube from the new anchor.

For same-face paths this is identical to `FastH3.gridPathCells` (single cube
call, no overhead).  For cross-face paths the greedy hops bridge the face
boundary with O(1) work per hop, then cube handles the remainder.
"""
function gridPathCellsRobust(start_::H3Index, end_::H3Index)
    err, cube_path = _cube_gridPathCells(start_, end_)
    if err == E_SUCCESS && !isempty(cube_path) && cube_path[end] == end_
        return cube_path
    end

    result = H3Index[]
    current = start_
    last_failed_bc = getBaseCellNumber(current)
    for _ in 1:10_000
        if current == end_
            if isempty(result) || result[end] != current
                push!(result, current)
            end
            return result
        end

        current_bc = getBaseCellNumber(current)
        if current_bc != last_failed_bc
            err, cube_path = _cube_gridPathCells(current, end_)
            if err == E_SUCCESS && !isempty(cube_path) && cube_path[end] == end_
                start_idx = (!isempty(result) && result[end] == cube_path[1]) ? 2 : 1
                append!(result, @view cube_path[start_idx:end])
                return result
            end
            last_failed_bc = current_bc
        end

        if isempty(result) || result[end] != current
            push!(result, current)
        end
        next = _step_toward(current, end_)
        next == current && break
        current = next
    end

    if isempty(result) || result[end] != end_
        push!(result, end_)
    end
    return result
end

"""
    gridDistanceRobust(origin::H3Index, h3::H3Index) -> (H3Error, Int64)

Number of grid hops along [`gridPathCellsRobust`](@ref): `length(path) - 1`.

For two valid cells at the same resolution, this matches [`FastH3.gridDistance`](@ref)
when the core implementation succeeds. When `gridDistance` returns `E_DOMAIN` or
other failures for valid same-res cells, a robust path is still computed.

Directed edges and other indices where `isValidCell` is false but
`gridDistance` succeeds are handled like `gridDistance` (no robust walk).
"""
function gridDistanceRobust(origin::H3Index, h3::H3Index)::Tuple{H3Error,Int64}
    err_gd, dist_gd = gridDistance(origin, h3)
    if err_gd == E_RES_MISMATCH || err_gd == E_CELL_INVALID
        return (err_gd, Int64(0))
    end
    if err_gd == E_SUCCESS && isValidCell(origin) && isValidCell(h3)
        path = gridPathCellsRobust(origin, h3)
        return (E_SUCCESS, Int64(length(path) - 1))
    end
    if err_gd == E_SUCCESS
        return (E_SUCCESS, dist_gd)
    end
    path = gridPathCellsRobust(origin, h3)
    return (E_SUCCESS, Int64(length(path) - 1))
end
