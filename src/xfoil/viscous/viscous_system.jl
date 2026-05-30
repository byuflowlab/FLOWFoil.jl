#=
Viscous-coupling influence matrices.

These matrices link the boundary-layer state (displacement thickness δ*
and edge velocity ue) to the inviscid γ distribution via the mass-injection
source representation. The end product `D` satisfies

    ue = ue_inv + D · (δ* · ue)

for the coupled boundary-layer / inviscid system. The pieces are:

- `B`      : Ψ at airfoil nodes from sources on airfoil and wake.
- `Bprime` : A⁻¹ B, sliced to drop the last (constant Ψ) row.
              Gives δγ at each airfoil node from each source.
- `Cgamma` : tangential wake-node velocity from airfoil γ.
- `Csigma` : tangential wake-node velocity from sources (airfoil and wake).
- `Dwake`  : tangential wake-node velocity from sources, = Cgamma·Bprime + Csigma.
- `Dprime` : d/ds matrix that converts (δ*·ue) on the surface to source.
- `D`      : `diag(d) · vcat(Bprime, Dwake) · Dprime`, the final operator
              used in the boundary-layer velocity-residual equation.

Sign conventions: FLOWFoil's A matrix is the negative of mfoil's, so the
B matrix is computed in mfoil convention here and then negated when
solving for Bprime (the negation cancels out, so Bprime is the same value
in both conventions). The downstream code consumes Bprime/Dwake/D in
mfoil's convention.
=#

"""
    create_B_influence_matrix(airfoil_panels, wake_panels)

Build the source-streamfunction influence matrix `B`. Each row `i`
(`1 ≤ i ≤ npanels + 1`) gives the streamfunction at airfoil node `i` from a
unit-strength source on:
- airfoil panel `j` (column `j`, `1 ≤ j ≤ npanels_airfoil`), constant source
- wake panel `k` (column `npanels_airfoil + k`, `1 ≤ k ≤ npanels_wake`),
  parameterized as piecewise-linear between wake-panel midpoints

The last row (Kutta condition) is left at zero. Output shape:
`(npanels+2) × (npanels_airfoil + npanels_wake)`.
"""
function create_B_influence_matrix(airfoil_panels, wake_panels)
    np_a = length(airfoil_panels)
    np_w = length(wake_panels)
    TF = promote_type(typeof(airfoil_panels[1].x1), typeof(wake_panels[1].x1))
    B = zeros(TF, np_a + 2, np_a + np_w)

    # x,y of each airfoil node (one more than panels)
    @inbounds for i in 1:(np_a + 1)
        if i ≤ np_a
            xi = airfoil_panels[i].x1
            yi = airfoil_panels[i].y1
        else
            xi = airfoil_panels[np_a].x2
            yi = airfoil_panels[np_a].y2
        end

        # ── airfoil-panel constant sources ───────────────────────────────────
        for j in 1:np_a
            xstar, ystar, _, _, logr1, logr2, β = mfoil_relative_geometry(xi, yi, airfoil_panels[j])
            B[i, j] = mfoil_constant_source_psi(
                xstar, ystar, logr1, logr2, airfoil_panels[j].length, β
            )
        end

        # ── wake-panel linear sources (piecewise linear between midpoints) ──
        for k in 1:np_w
            x1k = wake_panels[k].x1
            y1k = wake_panels[k].y1
            x2k = wake_panels[k].x2
            y2k = wake_panels[k].y2
            xmk = (x1k + x2k) / 2
            ymk = (y1k + y2k) / 2

            if k == np_w
                # Ghost extension on the last half-panel — σ continues constant
                x2g = 2 * x2k - xmk
                y2g = 2 * y2k - ymk
                l_l = sqrt((xmk - x1k)^2 + (ymk - y1k)^2)
                l_r = sqrt((x2g - xmk)^2 + (y2g - ymk)^2)
            else
                full = sqrt((x2k - x1k)^2 + (y2k - y1k)^2)
                l_l = full / 2
                l_r = full / 2
                x2g = x2k
                y2g = y2k
            end

            # left half-panel: from wake_node_k to wake_midpoint_k
            half_l = (;
                x1=x1k,
                y1=y1k,
                x2=xmk,
                y2=ymk,
                length=l_l,
                ctheta=(xmk - x1k) / l_l,
                stheta=(ymk - y1k) / l_l,
            )
            xstar_l, ystar_l, r1l, rml, logr1l, logrml, β_l = mfoil_relative_geometry(
                xi, yi, half_l
            )
            a_l, b_l = mfoil_linear_source_psi(
                xstar_l, ystar_l, logr1l, logrml, l_l, β_l, r1l, rml
            )

            # σ assignment for left half:
            #   k > 1 : left endpoint = (σ_{k-1} + σ_k)/2, right endpoint = σ_k
            #   k = 1 : left endpoint = 0 (TE start)        , right endpoint = σ_1
            if k > 1
                B[i, np_a + k] += a_l / 2 + b_l
                B[i, np_a + k - 1] += a_l / 2
            else
                B[i, np_a + k] += b_l
            end

            # right half-panel: from wake_midpoint_k to wake_node_{k+1}
            half_r = (;
                x1=xmk,
                y1=ymk,
                x2=x2g,
                y2=y2g,
                length=l_r,
                ctheta=(x2g - xmk) / l_r,
                stheta=(y2g - ymk) / l_r,
            )
            xstar_r, ystar_r, rmr, r2r, logrmr, logr2r, β_r = mfoil_relative_geometry(
                xi, yi, half_r
            )
            a_r, b_r = mfoil_linear_source_psi(
                xstar_r, ystar_r, logrmr, logr2r, l_r, β_r, rmr, r2r
            )

            # σ assignment for right half:
            #   k < npanels_wake : left = σ_k,         right = (σ_k + σ_{k+1})/2
            #   k = npanels_wake : left = σ_k,         right = σ_k (ghost)
            B[i, np_a + k] += a_r + b_r / 2
            if k < np_w
                B[i, np_a + k + 1] += b_r / 2
            else
                B[i, np_a + k] += b_r / 2
            end
        end
    end

    return B
end

"""
    create_Bprime_influence_matrix(B, A)

Solve the linear system `A · X = -B` (the sign flip absorbs the convention
mismatch between FLOWFoil's A and mfoil's B) and drop the last (constant
Ψ) row. Returns a `(npanels+1) × (npanels_airfoil + npanels_wake)` matrix.
"""
function create_Bprime_influence_matrix(B, A)
    Bprime = A \ (-B)
    
    return Bprime[1:(end - 1), :]
end

"""
    create_Cgamma_influence_matrix(airfoil_panels, te_geometry, wake_panels)

Tangential velocity at each wake node from each airfoil γ. Wake nodes are
the `npanels_wake + 1` panel boundaries; the tangent direction at each
node is taken from the wake panel that ends there (`tdir2` for the last
node, `tdir1` for all others). Returns shape
`(npanels_wake + 1) × (npanels_airfoil + 1)`.
"""
function create_Cgamma_influence_matrix(airfoil_panels, te_geometry, wake_panels)
    np_a = length(airfoil_panels)
    np_w = length(wake_panels)
    TF = promote_type(typeof(airfoil_panels[1].x1), typeof(wake_panels[1].x1))
    Cgamma = zeros(TF, np_w + 1, np_a + 1)

    has_blunt_te = te_geometry.length > 1e-12

    @inbounds for k in 1:(np_w + 1)
        if k ≤ np_w
            xk = wake_panels[k].x1
            yk = wake_panels[k].y1
            tdir = wake_panels[k].tdir1
        else
            xk = wake_panels[np_w].x2
            yk = wake_panels[np_w].y2
            tdir = wake_panels[np_w].tdir2
        end

        Vgx = zeros(TF, np_a + 1)
        Vgy = zeros(TF, np_a + 1)

        # airfoil linear-vortex panels
        for j in 1:np_a
            xstar, ystar, _, _, logr1, logr2, β = mfoil_relative_geometry(xk, yk, airfoil_panels[j])
            (ax, ay), (bx, by) = mfoil_linear_vortex_velocity(
                xstar, ystar, logr1, logr2, airfoil_panels[j].length, β, airfoil_panels[j]
            )
            Vgx[j] += ax
            Vgy[j] += ay
            Vgx[j + 1] += bx
            Vgy[j + 1] += by
        end

        # TE closure panel (source + vortex) — scaled by γ difference across the TE
        if has_blunt_te
            xstar, ystar, _, _, logr1, logr2, β = mfoil_relative_geometry(xk, yk, te_geometry)
            asx, asy = mfoil_constant_source_velocity(
                xstar, ystar, logr1, logr2, te_geometry.length, te_geometry
            )
            fs = te_geometry.tcp / 2
            Vgx[1] -= fs * asx
            Vgy[1] -= fs * asy
            Vgx[np_a + 1] += fs * asx
            Vgy[np_a + 1] += fs * asy

            (avx, avy), (bvx, bvy) = mfoil_linear_vortex_velocity(
                xstar, ystar, logr1, logr2, te_geometry.length, β, te_geometry
            )
            gx = (avx + bvx) * te_geometry.tdp / 2
            gy = (avy + bvy) * te_geometry.tdp / 2
            Vgx[1] += gx
            Vgy[1] += gy
            Vgx[np_a + 1] -= gx
            Vgy[np_a + 1] -= gy
        end

        # Project onto wake tangent
        for n in 1:(np_a + 1)
            Cgamma[k, n] = Vgx[n] * tdir[1] + Vgy[n] * tdir[2]
        end
    end

    return Cgamma
end

"""
    create_Csigma_influence_matrix(airfoil_panels, wake_panels)

Tangential velocity at each wake node from each source (airfoil panel
constant sources + wake panel linear sources). Returns shape
`(npanels_wake + 1) × (npanels_airfoil + npanels_wake)`.

Special cases:
- The first airfoil panel and the last airfoil panel do not contribute to
  the first wake node (those panels are adjacent to the TE and have their
  velocity treated implicitly).
- On the wake panel that contains the control point, an analytical
  self-influence formula is used in place of the panel integral.
- The ghost extension off the last wake panel is treated as a constant
  source.
"""
function create_Csigma_influence_matrix(airfoil_panels, wake_panels)
    np_a = length(airfoil_panels)
    np_w = length(wake_panels)
    TF = promote_type(typeof(airfoil_panels[1].x1), typeof(wake_panels[1].x1))
    Csigma = zeros(TF, np_w + 1, np_a + np_w)

    @inbounds for k in 1:(np_w + 1)
        if k ≤ np_w
            xk = wake_panels[k].x1
            yk = wake_panels[k].y1
            wake_tdir = wake_panels[k].tdir1
        else
            xk = wake_panels[np_w].x2
            yk = wake_panels[np_w].y2
            wake_tdir = wake_panels[np_w].tdir2
        end

        # airfoil constant sources
        for i in 1:np_a
            # First wake node is adjacent to first and last airfoil panels;
            # those panels' contributions are handled implicitly elsewhere.
            if k == 1 && (i == 1 || i == np_a)
                continue
            end
            p = airfoil_panels[i]
            xstar, ystar, _, _, logr1, logr2, β = mfoil_relative_geometry(xk, yk, p)
            vx, vy = mfoil_constant_source_velocity(xstar, ystar, logr1, logr2, p.length, p)
            Csigma[k, i] += vx * wake_tdir[1] + vy * wake_tdir[2]
        end

        # wake linear sources — piecewise linear between panel midpoints
        for j in 1:(np_w + 1)
            # half-panel endpoints around the j-th node
            if j == 1
                x1j = wake_panels[j].x1
                xmj = wake_panels[j].x1
                x2j = wake_panels[j].xbar
                y1j = wake_panels[j].y1
                ymj = wake_panels[j].y1
                y2j = wake_panels[j].ybar
            elseif j == np_w + 1
                x1j = wake_panels[j - 1].xbar
                xmj = wake_panels[j - 1].x2
                x2j = 2 * xmj - x1j
                y1j = wake_panels[j - 1].ybar
                ymj = wake_panels[j - 1].y2
                y2j = 2 * ymj - y1j
            else
                x1j = wake_panels[j - 1].xbar
                xmj = wake_panels[j].x1
                x2j = wake_panels[j].xbar
                y1j = wake_panels[j - 1].ybar
                ymj = wake_panels[j].y1
                y2j = wake_panels[j].ybar
            end

            l_l = sqrt((xmj - x1j)^2 + (ymj - y1j)^2)
            l_r = sqrt((x2j - xmj)^2 + (y2j - ymj)^2)

            half_l = (;
                x1=x1j,
                y1=y1j,
                x2=xmj,
                y2=ymj,
                length=l_l,
                ctheta=(xmj - x1j) / l_l,
                stheta=(ymj - y1j) / l_l,
            )
            half_r = (;
                x1=xmj,
                y1=ymj,
                x2=x2j,
                y2=y2j,
                length=l_r,
                ctheta=(x2j - xmj) / l_r,
                stheta=(y2j - ymj) / l_r,
            )

            if k == j
                # Self-influence — replaced with an analytical formula.
                if j == 1
                    l_lower = airfoil_panels[1].length
                    l_upper = airfoil_panels[end].length
                    Csigma[k, 1] += (log(l_lower / l_r) + 1) / (2π)
                    Csigma[k, np_a] += (log(l_upper / l_r) + 1) / (2π)
                    Csigma[k, np_a + 1] -= 1 / (2π)
                elseif j < np_w + 1
                    aa = log(l_l / l_r) / (4π)
                    Csigma[k, np_a + j - 1] += aa + 1 / (2π)
                    Csigma[k, np_a + j] += aa - 1 / (2π)
                end
                # j == np_w + 1 self-influence handled by the ghost extension below.
            else
                if j == 1
                    # only right half exists; left is degenerate
                    xstar_r, ystar_r, rmr, r2r, logrmr, logr2r, β_r = mfoil_relative_geometry(
                        xk, yk, half_r
                    )
                    (ax, ay), (bx, by) = mfoil_linear_source_velocity(
                        xstar_r, ystar_r, logrmr, logr2r, l_r, β_r, half_r
                    )
                    a_proj = ax * wake_tdir[1] + ay * wake_tdir[2]
                    b_proj = bx * wake_tdir[1] + by * wake_tdir[2]
                    Csigma[k, np_a + 1] += b_proj
                    Csigma[k, 1] += a_proj
                    Csigma[k, np_a] += a_proj
                elseif j == np_w + 1
                    # Ghost extension treated as a constant source over the
                    # full extended panel (x1j, y1j) → (x2j, y2j).
                    l_end = sqrt((x2j - x1j)^2 + (y2j - y1j)^2)
                    panel_end = (;
                        x1=x1j,
                        y1=y1j,
                        x2=x2j,
                        y2=y2j,
                        length=l_end,
                        ctheta=(x2j - x1j) / l_end,
                        stheta=(y2j - y1j) / l_end,
                    )
                    xstar_e, ystar_e, _, _, logr1e, logr2e, _ = mfoil_relative_geometry(
                        xk, yk, panel_end
                    )
                    vx, vy = mfoil_constant_source_velocity(
                        xstar_e, ystar_e, logr1e, logr2e, l_end, panel_end
                    )
                    Csigma[k, np_a + np_w] += vx * wake_tdir[1] + vy * wake_tdir[2]
                else
                    xstar_l, ystar_l, r1l, rml, logr1l, logrml, β_l = mfoil_relative_geometry(
                        xk, yk, half_l
                    )
                    xstar_r, ystar_r, rmr, r2r, logrmr, logr2r, β_r = mfoil_relative_geometry(
                        xk, yk, half_r
                    )
                    (a1x, a1y), (b1x, b1y) = mfoil_linear_source_velocity(
                        xstar_l, ystar_l, logr1l, logrml, l_l, β_l, half_l
                    )
                    (a2x, a2y), (b2x, b2y) = mfoil_linear_source_velocity(
                        xstar_r, ystar_r, logrmr, logr2r, l_r, β_r, half_r
                    )

                    a1p = a1x * wake_tdir[1] + a1y * wake_tdir[2]
                    b1p = b1x * wake_tdir[1] + b1y * wake_tdir[2]
                    a2p = a2x * wake_tdir[1] + a2y * wake_tdir[2]
                    b2p = b2x * wake_tdir[1] + b2y * wake_tdir[2]
                    Csigma[k, np_a + j - 1] += a1p + b1p / 2
                    Csigma[k, np_a + j] += a2p / 2 + b2p
                end
            end
        end
    end

    return Csigma
end

"""
    create_Dwake_influence_matrix(Bprime, Cgamma, Csigma)

Tangential wake-node velocity from sources, combining the γ-mediated
contribution (`Cgamma·Bprime`) with the direct source contribution
(`Csigma`). The first row is overwritten with the last row of `Bprime`
since the first wake node coincides with the TE.
"""
function create_Dwake_influence_matrix(Bprime, Cgamma, Csigma)
    Dwake = Cgamma * Bprime + Csigma
    Dwake[1, :] = Bprime[end, :]
    
    return Dwake
end

"""
    create_Dprime_influence_matrix(s_airfoil, s_wake, d_airfoil, d=nothing)

Sparse `(npanels_airfoil + npanels_wake) × (npanels_airfoil + npanels_wake + 2)`
matrix that converts a node-wise `(δ*·ue)` distribution into per-panel
mass-flux sources via finite-difference d/ds. The directional sign
vector `d_airfoil` carries the airfoil γ sign at each node so that the
wake "downstream" direction is unambiguous.
"""
function create_Dprime_influence_matrix(s_airfoil, s_wake, d_airfoil, d=nothing)
    np_a = length(s_airfoil) - 1
    np_w = length(s_wake) - 1
    Dprime = spzeros(np_a + np_w, np_a + np_w + 2)

    for i in 1:np_a
        ds = s_airfoil[i + 1] - s_airfoil[i]
        Dprime[i, i] = -d_airfoil[i] / ds
        Dprime[i, i + 1] = d_airfoil[i + 1] / ds
    end

    for k in 1:np_w
        ds = s_wake[k + 1] - s_wake[k]
        Dprime[np_a + k, np_a + 1 + k] = -1 / ds
        Dprime[np_a + k, np_a + 1 + k + 1] = 1 / ds
    end

    return Dprime
end

"""
    get_d_directional_vectors(numpanels_wake, gamma_at_alpha)

Per-node sign vector for the airfoil + wake. On the airfoil it is the
sign of the tangential velocity (so that the BL solves work in the
positive-flow direction on each surface); on the wake it is +1 everywhere
(wake direction is always downstream of the TE).
"""
function get_d_directional_vectors(numpanels_wake, gamma_at_alpha)
    d_airfoil = sign.(gamma_at_alpha)
    d_wake = ones(eltype(gamma_at_alpha), numpanels_wake + 1)
    d = vcat(d_airfoil, d_wake)
    
    return d_airfoil, d_wake, d
end

"""
    create_D_influence_matrix(airfoil_panels, te_geometry, wake_panels, A, gamma_at_alpha)

Build all the viscous-coupling influence matrices and the final operator
`D = diag(d) · vcat(Bprime, Dwake) · Dprime`. Returns `(D, d, Bprime, Dwake)`.

# Arguments
- `airfoil_panels` : Vector of per-panel NamedTuples (`build_viscous_airfoil_panels`).
- `te_geometry`    : TE closure-panel NamedTuple (`build_te_geometry`).
- `wake_panels`    : Vector of wake panel NamedTuples (`build_wake`).
- `A`              : `(npanels+2)×(npanels+2)` inviscid influence matrix from FLOWFoil.
- `gamma_at_alpha` : Surface tangential γ at the current angle of attack
                     (vector of length `npanels + 1`).
"""
function create_D_influence_matrix(
    airfoil_panels, te_geometry, wake_panels, A, gamma_at_alpha
)
    B = create_B_influence_matrix(airfoil_panels, wake_panels)
    Bprime = create_Bprime_influence_matrix(B, A)
    Cgamma = create_Cgamma_influence_matrix(airfoil_panels, te_geometry, wake_panels)
    Csigma = create_Csigma_influence_matrix(airfoil_panels, wake_panels)
    Dwake = create_Dwake_influence_matrix(Bprime, Cgamma, Csigma)

    np_w = length(wake_panels)
    d_airfoil, _, d = get_d_directional_vectors(np_w, gamma_at_alpha)

    s_airfoil = vcat(0.0, cumsum([p.length for p in airfoil_panels]))
    s_wake = vcat(s_airfoil[end], s_airfoil[end] .+ cumsum([p.length for p in wake_panels]))
    Dprime = create_Dprime_influence_matrix(s_airfoil, s_wake, d_airfoil)

    D = spdiagm(d) * vcat(Bprime, Dwake) * Dprime
    
    return D, d, Bprime, Dwake
end
