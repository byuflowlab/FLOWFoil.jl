#=
Wake construction for the viscous coupling.

Builds the streamline wake by integrating the inviscid velocity field
starting just behind the trailing edge. Each step is a predictor-corrector
in the local velocity direction. Returns a Vector of NamedTuples (one per
wake panel) plus a `(numpanels_wake+1) × 2` matrix of reference wake-edge
velocities (one column per unit cos α / sin α component) and the
wake-gap profile that continues the TE thickness downstream.

Also provides `induced_velocity_at_point`, the inviscid velocity field
evaluated at an arbitrary point, which the wake march and other parts of
the coupling code use.
=#

"""
    induced_velocity_at_point(airfoil_panels, te_geometry, gammas, point, Vinf, cosalpha, sinalpha)

Inviscid velocity at `point = (x, y)` produced by:
- linear vortex panels along the airfoil with strengths `gammas` (length
  `npanels + 1`, one per node),
- a constant-source + linear-vortex contribution from the TE closure
  panel scaled by `(γ_end − γ_1)·tcp/2` and `(γ_1 − γ_end)·tdp/2`,
- plus the freestream `Vinf · (cosalpha, sinalpha)`.
"""
function induced_velocity_at_point(
    airfoil_panels, te_geometry, gammas, point, Vinf, cosalpha, sinalpha
)
    Vx = zero(eltype(gammas))
    Vy = zero(eltype(gammas))
    xi, yi = point

    @inbounds for j in 1:length(airfoil_panels)
        p = airfoil_panels[j]
        xstar, ystar, _, _, logr1, logr2, β = mfoil_relative_geometry(xi, yi, p)
        (ax, ay), (bx, by) = mfoil_linear_vortex_velocity(
            xstar, ystar, logr1, logr2, p.length, β, p
        )
        Vx += ax * gammas[j] + bx * gammas[j + 1]
        Vy += ay * gammas[j] + by * gammas[j + 1]
    end

    # TE closure-panel source-vortex contributions (sharp TE: skip).
    if te_geometry.length > 0
        xstar, ystar, _, _, logr1, logr2, β = mfoil_relative_geometry(xi, yi, te_geometry)
        asx, asy = mfoil_constant_source_velocity(
            xstar, ystar, logr1, logr2, te_geometry.length, te_geometry
        )
        factor_s = te_geometry.tcp * (gammas[end] - gammas[1]) / 2
        Vx += factor_s * asx
        Vy += factor_s * asy

        (avx, avy), (bvx, bvy) = mfoil_linear_vortex_velocity(
            xstar, ystar, logr1, logr2, te_geometry.length, β, te_geometry
        )
        factor_v = te_geometry.tdp * (gammas[1] - gammas[end]) / 2
        Vx += factor_v * (avx + bvx)
        Vy += factor_v * (avy + bvy)
    end

    # Freestream
    Vx += Vinf * cosalpha
    Vy += Vinf * sinalpha

    return Vx, Vy
end

# ─── wake build ──────────────────────────────────────────────────────────────
"""
    build_wake(airfoil_panels, te_geometry, s_airfoil, gamma_ref, alpha,
               Vinf, chord, wake_length, numpanels_wake)

Construct the airfoil wake by integrating along the inviscid streamline
starting just behind the TE midpoint. Returns:

- `wake_panels` — Vector of NamedTuples with fields `x1, y1, x2, y2, xbar,
  ybar, length, stheta, ctheta, tdir1, tdir2` (last two are the unit
  velocity tangent at each panel endpoint, used by the BL initialization
  and wake gap continuation).
- `ue_wake_ref` — `(numpanels_wake+1) × 2` matrix; column 1 is the
  tangential wake-edge velocity for unit cos α (with sin α = 0), column 2
  is for unit sin α (cos α = 0).
- `wake_gap` — `(numpanels_wake+1)`-vector; TE-gap continuation profile
  used by the BL deltastar residuals on the wake.

# Arguments
- `airfoil_panels` : Vector of airfoil panel NamedTuples (see `build_viscous_airfoil_panels`).
- `te_geometry`    : TE closure-panel NamedTuple (see `build_te_geometry`).
- `s_airfoil`      : arc-length coordinate along the airfoil (`length npanels + 1`).
- `gamma_ref`      : `(npanels + 1) × 2` matrix of node vortex strengths
                     for unit cos α (col 1) and unit sin α (col 2).
- `alpha`          : angle of attack in **radians**.
- `Vinf`, `chord`  : freestream speed and reference chord.
- `wake_length`    : wake length in chords.
- `numpanels_wake` : number of wake panels.
"""
function build_wake(
    airfoil_panels,
    te_geometry,
    s_airfoil,
    gamma_ref,
    alpha,
    Vinf,
    chord,
    wake_length,
    numpanels_wake,
)
    cosalpha = cos(alpha)
    sinalpha = sin(alpha)
    γs = gamma_ref * [cosalpha, sinalpha]

    # Wake node spacing — first interval matches local airfoil panel scale near TE
    ds1 = max((s_airfoil[2] - s_airfoil[1] + s_airfoil[end] - s_airfoil[end - 1]) / 2, 1e-3)
    sv = geometric_spacing(ds1, wake_length * chord, numpanels_wake + 1)

    # First wake node: just behind the TE midpoint along the TE-gap normal
    xy1 = (airfoil_panels[1].x1, airfoil_panels[1].y1)
    xyN = (airfoil_panels[end].x2, airfoil_panels[end].y2)
    nx_TE = xyN[1] - xy1[1]
    ny_TE = xyN[2] - xy1[2]
    tx_TE = ny_TE       # 90° rotation: tangent vector
    ty_TE = -nx_TE
    TE_mid_x = (xy1[1] + xyN[1]) / 2
    TE_mid_y = (xy1[2] + xyN[2]) / 2

    TF = promote_type(eltype(s_airfoil), eltype(gamma_ref))
    wake_nodes_x = zeros(TF, numpanels_wake + 1)
    wake_nodes_y = zeros(TF, numpanels_wake + 1)
    tdir_x = zeros(TF, numpanels_wake + 1)
    tdir_y = zeros(TF, numpanels_wake + 1)

    wake_nodes_x[1] = TE_mid_x + 1e-5 * tx_TE * chord
    wake_nodes_y[1] = TE_mid_y + 1e-5 * ty_TE * chord

    # Predictor-corrector streamline integration
    for i in 1:numpanels_wake
        v1x, v1y = induced_velocity_at_point(
            airfoil_panels,
            te_geometry,
            γs,
            (wake_nodes_x[i], wake_nodes_y[i]),
            Vinf,
            cosalpha,
            sinalpha,
        )
        m1 = sqrt(v1x^2 + v1y^2)
        v1x /= m1
        v1y /= m1
        tdir_x[i] = v1x
        tdir_y[i] = v1y

        ds = sv[i + 1] - sv[i]
        # predictor
        px = wake_nodes_x[i] + ds * v1x
        py = wake_nodes_y[i] + ds * v1y

        v2x, v2y = induced_velocity_at_point(
            airfoil_panels, te_geometry, γs, (px, py), Vinf, cosalpha, sinalpha
        )
        m2 = sqrt(v2x^2 + v2y^2)
        v2x /= m2
        v2y /= m2
        tdir_x[i + 1] = v2x
        tdir_y[i + 1] = v2y

        # corrector — average of predictor/corrector directions
        wake_nodes_x[i + 1] = wake_nodes_x[i] + ds * (v1x + v2x) / 2
        wake_nodes_y[i + 1] = wake_nodes_y[i] + ds * (v1y + v2y) / 2
    end

    # Reference wake-edge velocities (unit cos α and unit sin α components)
    ue_wake_ref = zeros(TF, numpanels_wake + 1, 2)
    for i in 1:(numpanels_wake + 1)
        v0x, v0y = induced_velocity_at_point(
            airfoil_panels,
            te_geometry,
            gamma_ref[:, 1],
            (wake_nodes_x[i], wake_nodes_y[i]),
            Vinf,
            1.0,
            0.0,
        )
        v90x, v90y = induced_velocity_at_point(
            airfoil_panels,
            te_geometry,
            gamma_ref[:, 2],
            (wake_nodes_x[i], wake_nodes_y[i]),
            Vinf,
            0.0,
            1.0,
        )
        ue_wake_ref[i, 1] = v0x * tdir_x[i] + v0y * tdir_y[i]
        ue_wake_ref[i, 2] = v90x * tdir_x[i] + v90y * tdir_y[i]
    end

    # Assemble wake panel NamedTuples
    wake_panels = [
        let
            x1 = wake_nodes_x[i]
            y1 = wake_nodes_y[i]
            x2 = wake_nodes_x[i + 1]
            y2 = wake_nodes_y[i + 1]
            len = sqrt((x2 - x1)^2 + (y2 - y1)^2)
            (;
                x1=x1,
                y1=y1,
                x2=x2,
                y2=y2,
                xbar=(x1 + x2) / 2,
                ybar=(y1 + y2) / 2,
                length=len,
                stheta=(y2 - y1) / len,
                ctheta=(x2 - x1) / len,
                tdir1=(tdir_x[i], tdir_y[i]),
                tdir2=(tdir_x[i + 1], tdir_y[i + 1]),
            )
        end for i in 1:numpanels_wake
    ]

    # Wake-gap continuation profile
    wake_gap = zeros(TF, numpanels_wake + 1)
    panel_lengths = [p.length for p in wake_panels]
    f_length = 2.5
    dtdx = clamp(te_geometry.dtdx, -3 / f_length, 3 / f_length)
    Lw = f_length * te_geometry.hTE
    wake_gap[1] = te_geometry.hTE
    s_along_wake = cumsum(panel_lengths)
    for i in 2:(numpanels_wake + 1)
        ξ = s_along_wake[i - 1] / Lw
        if ξ < 1
            wake_gap[i] = te_geometry.hTE * (1 + (2 + f_length * dtdx) * ξ) * (1 - ξ)^2
        end
    end

    return wake_panels, ue_wake_ref, wake_gap
end
