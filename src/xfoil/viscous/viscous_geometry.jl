#=
Geometry helpers used by the viscous coupling. 
=#

"""
    build_viscous_airfoil_panels(panel_geometry)

Build a `Vector` of per-panel NamedTuples for the airfoil from FLOWFoil's
single-body `panel_geometry`. Each entry contains: `x1, y1, x2, y2` (panel
endpoints), `xbar, ybar` (midpoint), `length` (panel length), and
`stheta, ctheta` (panel tangent direction sin/cos).

Assumes nodes are ordered TE → lower → LE → upper → TE (Hess-Smith / xfoil
convention).
"""
function build_viscous_airfoil_panels(panel_geometry)
    npanels = panel_geometry.npanels
    
    return [
        let
            x1 = panel_geometry.panel_edges[i, 1, 1]
            y1 = panel_geometry.panel_edges[i, 1, 2]
            x2 = panel_geometry.panel_edges[i, 2, 1]
            y2 = panel_geometry.panel_edges[i, 2, 2]
            len = panel_geometry.panel_lengths[i]
            stheta = (y2 - y1) / len
            ctheta = (x2 - x1) / len
            (;
                x1=x1,
                y1=y1,
                x2=x2,
                y2=y2,
                xbar=(x1 + x2) / 2,
                ybar=(y1 + y2) / 2,
                length=len,
                stheta=stheta,
                ctheta=ctheta,
            )
        end for i in 1:npanels
    ]
end

"""
    build_te_geometry(airfoil_panels)

Build the trailing-edge gap panel descriptor used by the wake and influence
matrices. Returns a NamedTuple with the same endpoint/midpoint/length/direction
fields as a regular panel plus the TE-specific scalars `tdp`, `tcp`, `hTE`, `dtdx`.

# TE scalars
- `tdp` : dot product of the TE bisection unit vector with the TE gap vector.
- `tcp` : 2D cross of the same two vectors.
- `hTE` : signed TE gap length normal to the bisector (used to set up the wake gap).
- `dtdx` : signed change in panel direction across the TE (TE thickness slope).
"""
function build_te_geometry(airfoil_panels)
    x1 = airfoil_panels[end].x2
    y1 = airfoil_panels[end].y2
    x2 = airfoil_panels[1].x1
    y2 = airfoil_panels[1].y1

    len = sqrt((x2 - x1)^2 + (y2 - y1)^2)

    # Sharp TE: TE_upper ≈ TE_lower. The TE source-vortex closure panel
    # contributes nothing, so return zeros for the gap-dependent fields.
    if len < 1e-12
        TF = typeof(x1)
        return (;
            x1=x1, y1=y1, x2=x2, y2=y2,
            xbar=(x1 + x2) / 2, ybar=(y1 + y2) / 2,
            length=zero(TF),
            stheta=zero(TF),
            ctheta=zero(TF),
            tdp=zero(TF),
            tcp=zero(TF),
            hTE=zero(TF),
            dtdx=zero(TF),
        )
    end

    stheta = (y2 - y1) / len
    ctheta = (x2 - x1) / len

    # Tangent direction at each side of the TE
    t̂TE, t1, t2 = _te_tangent_bisector(
        x2,
        y2,
        airfoil_panels[2].x1,
        airfoil_panels[2].y1,
        airfoil_panels[end].x1,
        airfoil_panels[end].y1,
        x1,
        y1,
    )
    p̂TE = _te_gap_unit_vector(x2, y2, x1, y1)
    tdp = t̂TE[1] * p̂TE[1] + t̂TE[2] * p̂TE[2]
    tcp = t̂TE[1] * p̂TE[2] - t̂TE[2] * p̂TE[1]

    sx = x1 - x2
    sy = y1 - y2
    hTE = -sx * t̂TE[2] + sy * t̂TE[1]
    dtdx = t1[1] * t2[2] - t2[1] * t1[2]

    return (;
        x1=x1,
        y1=y1,
        x2=x2,
        y2=y2,
        xbar=(x1 + x2) / 2,
        ybar=(y1 + y2) / 2,
        length=len,
        stheta=stheta,
        ctheta=ctheta,
        tdp=tdp,
        tcp=tcp,
        hTE=hTE,
        dtdx=dtdx,
    )
end

function _te_tangent_bisector(x1, y1, x2, y2, xn, yn, xnp1, ynp1)
    v1x = x1 - x2
    v1y = y1 - y2
    n1 = sqrt(v1x^2 + v1y^2)
    v1x /= n1
    v1y /= n1

    vnx = xnp1 - xn
    vny = ynp1 - yn
    nn = sqrt(vnx^2 + vny^2)
    vnx /= nn
    vny /= nn

    bx = v1x + vnx
    by = v1y + vny
    nb = sqrt(bx^2 + by^2)

    return (bx / nb, by / nb), (v1x, v1y), (vnx, vny)
end

function _te_gap_unit_vector(x1, y1, xnp1, ynp1)
    vx = xnp1 - x1
    vy = ynp1 - y1
    n = sqrt(vx^2 + vy^2)
    return (vx / n, vy / n)
end

"""
    find_stagnation_split(x, y, Vt, airfoil_panels, xft_l, xft_u)

Find the stagnation point on the airfoil from the tangential-velocity
distribution `Vt` (one value per node) and split the arc-length coordinate
into upper and lower surfaces.

# Arguments
- `x, y` : node coordinates (TE → lower → LE → upper → TE).
- `Vt`   : signed surface tangential velocity at each node (Hess-Smith
           convention: negative on lower surface, positive on upper).
- `airfoil_panels` : Vector of panel NamedTuples from `build_viscous_airfoil_panels`.
- `xft_l, xft_u`   : forced-transition x/c on lower/upper surfaces. If the
                     value is ≥ 1, no forcing is applied on that surface.

# Returns
NamedTuple with:
- `s`       : arc-length coordinate at each node, starting at TE.
- `su`      : arc-length on the upper surface measured from the stagnation point.
- `idx_u`   : `[i_stagnation, i_TE_upper]` node indices for the upper surface.
- `sl`      : arc-length on the lower surface measured from the stagnation point (positive going downstream toward TE).
- `idx_l`   : `[i_stagnation - 1, 1]` node indices for the lower surface.
- `s_stag`  : arc-length location of the stagnation point.
- `sft_l`   : forced-transition arc-length on lower surface (set to `s[end]` when not forced).
- `sft_u`   : forced-transition arc-length on upper surface (set to `s[end]` when not forced).
"""
function find_stagnation_split(x, y, Vt, airfoil_panels, xft_l, xft_u)
    idx = findfirst(z -> z > 0, Vt)
    if idx === nothing || idx == 1
        error("Could not locate stagnation point in Vt distribution (no sign change found).")
    end

    idx_l = [idx - 1, 1]
    idx_u = [idx, length(Vt)]

    lengths = [p.length for p in airfoil_panels]
    s = vcat(0.0, cumsum(lengths))

    s_stag = ImplicitAD.implicit(
        _solve_stagnation_state_implicit, _residual_stagnation_implicit, vcat(s, Vt), nothing
    )

    sl = -(s[(idx - 1):-1:1] .- s_stag)
    su = s[idx:end] .- s_stag

    sft_l = xft_l < 1 ? FLOWMath.linear(x[(idx - 1):-1:1], sl, xft_l) : s[end]
    sft_u = xft_u < 1 ? FLOWMath.linear(x[idx:end], su, xft_u) : s[end]

    return (;
        s=s,
        su=su,
        idx_u=idx_u,
        sl=sl,
        idx_l=idx_l,
        s_stag=s_stag,
        sft_l=sft_l,
        sft_u=sft_u,
    )
end

function _solve_stagnation_state_implicit(xstate, _)
    npts = length(xstate) ÷ 2
    s = xstate[1:npts]
    Vt = xstate[(npts + 1):end]

    spline = FLOWMath.Akima(s, Vt)
    s_guess = FLOWMath.akima(Vt, s, 0.0)

    return Roots.find_zero(spline, s_guess)
end

function _residual_stagnation_implicit(y, xstate, _)
    npts = length(xstate) ÷ 2
    s = xstate[1:npts]
    Vt = xstate[(npts + 1):end]

    spline = FLOWMath.Akima(s, Vt)
    
    return spline(y)
end

"""
    resplit_at_stagnation(s, ue, idx_u, idx_l, idx_w, swake0)

Re-split the arc-length coordinate into upper/lower/wake segments around an
updated stagnation point implied by the current `ue` state vector. The
stagnation point is located by linear interpolation of the signed
`ue` distribution to zero.

# Returns
- `su, sl, sw` : per-surface arc-length vectors measured from the new stagnation point.
- `s_stag`     : new stagnation point arc-length.
"""
function resplit_at_stagnation(s, ue, idx_u, idx_l, idx_w, swake0)
    Vt = vcat(reverse(-ue[idx_l[1]:-1:idx_l[end]]), ue[idx_u[1]:idx_u[end]])

    s_stag = FLOWMath.linear(Vt, s, 0.0)
    su = s[idx_u[1]:idx_u[end]] .- s_stag
    sl = -(s[idx_l[1]:-1:idx_l[end]] .- s_stag)
    sw = swake0 .- s_stag

    return su, sl, sw, s_stag
end

"""
    arclength_vectors(airfoil_panels, wake_panels, s_stag=0.0)

Build the arc-length coordinate along the airfoil (`s_airfoil`) and along
the wake measured from the stagnation point (`s_wake`), plus a copy
`swake0` of `s_wake` referenced to the airfoil start.
"""
function arclength_vectors(airfoil_panels, wake_panels, s_stag=0.0)
    s_airfoil = vcat(0.0, cumsum([p.length for p in airfoil_panels]))
    swake0 = vcat(s_airfoil[end], s_airfoil[end] .+ cumsum([p.length for p in wake_panels]))
    s_wake = swake0 .- s_stag
    
    return s_airfoil, s_wake, swake0
end

"""
    geometric_spacing(dx0, L, np)

Return an `np`-length vector of geometrically-spaced points on `[0, L]`
with first interval `dx0`. Used for wake spacing.
"""
function geometric_spacing(dx0, L, np)
    N = np - 1
    d = L / dx0

    a = N * (N - 1) * (N - 2) / 6
    b = N * (N - 1) / 2
    c = N - d

    disc = max(b^2 - 4a * c, 0.0)
    r = 1 + (-b + sqrt(disc)) / (2a)

    for _ in 1:10
        R = r^N - 1 - d * (r - 1)
        Rr = N * r^(N - 1) - d
        dr = -R / Rr
        r += dr
        if abs(dr) < 1e-6
            break
        end
    end

    dx = dx0 .* r .^ (0:(N - 1))
    
    return vcat(0.0, cumsum(dx))
end
