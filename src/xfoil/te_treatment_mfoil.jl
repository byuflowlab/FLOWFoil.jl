#=
Blunt-trailing-edge A-matrix contribution, implemented in the mfoil/xfoil
convention. The existing FLOWFoil treatment in `assemble_influence_matrix`
gave cl values that disagree with mfoil and the legacy Fortran Xfoil on
blunt-TE airfoils (e.g. NACA 2412 at α=5° gives cl≈0.77 vs the reference
≈0.858). This module supplies a self-contained TE contribution that
matches mfoil exactly while preserving FLOWFoil's overall A sign
convention (which is the negative of mfoil's — both b and A are sign-flipped
between the two libraries, so γ is unchanged for sharp-TE airfoils).

The math: at airfoil node `i`, the blunt-TE source-vortex closure panel
contributes a streamfunction influence

    Ψ_TE(i) = (γ_lower - γ_upper)/2 · (|tdp|·Ψ_γ_panel + |tcp|·Ψ_σ_panel)

where Ψ_γ_panel is the constant-vortex streamfunction integral (I1 in
mfoil's notation) over the TE gap panel and Ψ_σ_panel is the constant-
source streamfunction integral (I3 with branch-cut shift). `tdp` is the
dot product of the TE bisector unit vector and the TE gap unit vector;
`tcp` is the corresponding (signed) 2D cross product. The TE panel is
oriented in mfoil's convention: node 1 = TE_upper, node 2 = TE_lower.
=#

"""
    apply_blunt_te_treatment_mfoil!(amat, nidx, nodes)

Add the mfoil-style blunt trailing-edge source-vortex contribution to the
influence matrix `amat` for the body whose node range is `nidx` and whose
ordered node coordinates are the rows `nodes[nidx, :]`. Uses FLOWFoil's
overall A sign convention (which is negated from mfoil's), so the
contribution applied here is the *negative* of mfoil's direct formula.
"""
function apply_blunt_te_treatment_mfoil!(amat, nidx, nodes)
    # TE panel endpoints in mfoil convention: node1 = TE_upper, node2 = TE_lower
    i_lower = nidx[1]
    i_upper = nidx[end]
    xu = nodes[i_upper, 1]
    yu = nodes[i_upper, 2]
    xl = nodes[i_lower, 1]
    yl = nodes[i_lower, 2]

    lj = sqrt((xl - xu)^2 + (yl - yu)^2)
    if lj == 0
        return amat   # sharp TE; should not be reached, but guard anyway
    end
    ctheta = (xl - xu) / lj
    stheta = (yl - yu) / lj

    # TE bisector and gap unit vectors
    # bisector points AFT (away from LE)
    x_first_next = nodes[i_lower + 1, 1]
    y_first_next = nodes[i_lower + 1, 2]
    x_last_prev = nodes[i_upper - 1, 1]
    y_last_prev = nodes[i_upper - 1, 2]

    v1x = xl - x_first_next
    v1y = yl - y_first_next
    n1 = sqrt(v1x^2 + v1y^2)
    v1x /= n1
    v1y /= n1

    vnx = xu - x_last_prev
    vny = yu - y_last_prev
    nn = sqrt(vnx^2 + vny^2)
    vnx /= nn
    vny /= nn

    bx = v1x + vnx
    by = v1y + vny
    nb = sqrt(bx^2 + by^2)
    bxh = bx / nb
    byh = by / nb

    # gap vector from TE_lower to TE_upper (mfoil sign convention)
    gx = (xu - xl) / lj
    gy = (yu - yl) / lj

    tdp = bxh * gx + byh * gy
    tcp = bxh * gy - byh * gx

    abs_tdp = abs(tdp)
    abs_tcp = abs(tcp)

    # Loop over field nodes on this body
    for i in nidx
        xi = nodes[i, 1]
        yi = nodes[i, 2]

        # Panel-local coords of node i (TE panel goes TE_upper -> TE_lower)
        dx = xi - xu
        dy = yi - yu
        xstar = dx * ctheta + dy * stheta
        ystar = -dx * stheta + dy * ctheta

        # Distances from node i to TE_upper and TE_lower
        r1 = sqrt(dx^2 + dy^2)                       # distance to TE_upper (panel node 1)
        r2 = sqrt((xi - xl)^2 + (yi - yl)^2)         # distance to TE_lower (panel node 2)
        logr1 = r1 == 0 ? 0.0 : log(r1)
        logr2 = r2 == 0 ? 0.0 : log(r2)

        # Subtended angle at node i, with vertex at the panel
        β = _te_beta(xi, yi, xu, yu, xl, yl)

        # Linear-vortex psi integral I1 (mfoil convention)
        I1 = (xstar * (logr1 - logr2) + lj * logr2 + ystar * β - lj) / (2π)

        # Constant-source psi integral (mfoil convention) with branch-cut shift
        θ1 = atan(ystar, xstar)
        θ2 = atan(ystar, xstar - lj)
        if r1 == 0
            θ1 = π
            θ2 = π
            β_eff = 0.0
        elseif r2 == 0
            θ1 = 0.0
            θ2 = 0.0
            β_eff = 0.0
        else
            β_eff = β
        end
        I3 = (ystar * (logr1 - logr2) - xstar * β_eff + lj * θ2) / (2π)
        Ψσ = (θ1 + θ2) > π ? (I3 - 0.25 * lj) : (I3 + 0.75 * lj)

        # mfoil applies: A[i, lower] += +0.5*(|tdp|*I1 + |tcp|*Ψσ)
        #                A[i, upper] += -0.5*(|tdp|*I1 + |tcp|*Ψσ)
        # FLOWFoil's overall A is sign-flipped from mfoil's, so we apply the negative.
        contribution = 0.5 * (abs_tdp * I1 + abs_tcp * Ψσ)
        amat[i, nidx[1]]   -= contribution
        amat[i, nidx[end]] += contribution
    end

    return amat
end

# Signed angle subtended at field point (xi, yi) by the line segment from
# (x1, y1) to (x2, y2). Matches mfoil's `get_beta_mfoil`.
function _te_beta(xi, yi, x1, y1, x2, y2)
    if (xi == x1 && yi == y1) || (xi == x2 && yi == y2)
        return 0.0
    end
    num = (x1 - xi) * (y2 - yi) - (y1 - yi) * (x2 - xi)
    den = (x1 - xi) * (x2 - xi) + (y1 - yi) * (y2 - yi)
    
    return atan(num, den)
end
