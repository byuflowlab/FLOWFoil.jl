#=
mfoil-based panel-integral helpers (streamfunction and velocity influences)
used by the viscous coupling: B / Bprime / Cgamma / Csigma / Dwake / D matrices,
the wake march, and the boundary-layer residual setup.

A "panel" here is any NamedTuple with fields
`x1, y1, x2, y2, ctheta, stheta, length`.
=#

# Signed angle subtended at field point (xi, yi) by the segment (x1,y1)→(x2,y2).
function mfoil_beta(xi, yi, x1, y1, x2, y2)
    if (xi == x1 && yi == y1) || (xi == x2 && yi == y2)
        return zero(promote_type(typeof(xi), typeof(yi), typeof(x1)))
    end
    num = (x1 - xi) * (y2 - yi) - (y1 - yi) * (x2 - xi)
    den = (x1 - xi) * (x2 - xi) + (y1 - yi) * (y2 - yi)
    
    return atan(num, den)
end

# Panel-local coordinates of the field point (xi, yi), log-distances to the
# panel endpoints, and the angle subtended at the field point.
function mfoil_relative_geometry(xi, yi, panel)
    dx1 = xi - panel.x1
    dy1 = yi - panel.y1
    dx2 = xi - panel.x2
    dy2 = yi - panel.y2

    xstar = dx1 * panel.ctheta + dy1 * panel.stheta
    ystar = -dx1 * panel.stheta + dy1 * panel.ctheta

    r1 = sqrt(dx1^2 + dy1^2)
    r2 = sqrt(dx2^2 + dy2^2)
    logr1 = r1 == 0 ? zero(r1) : log(r1)
    logr2 = r2 == 0 ? zero(r2) : log(r2)
    β = (r1 == 0 || r2 == 0) ? zero(xi) : mfoil_beta(xi, yi, panel.x1, panel.y1, panel.x2, panel.y2)

    return xstar, ystar, r1, r2, logr1, logr2, β
end

# Linear-vortex streamfunction influence integrals over a panel, evaluated at a
# field point. Returns (I1, I2) such that the contribution to Ψ from a panel
# with γ varying linearly from γ_left to γ_right is `I1·γ_left + I2·γ_right`.
function mfoil_linear_vortex_psi(xstar, ystar, logr1, logr2, lj, β, r1, r2)
    I1 = (xstar * (logr1 - logr2) + lj * logr2 + ystar * β - lj) / (2π)
    I2 = (xstar * I1 + (r2^2 * (logr2 - 0.5) - r1^2 * (logr1 - 0.5)) / (4π)) / lj
    # `a` multiplies γ_left, `b` multiplies γ_right
    
    return I1 - I2, I2
end

# Constant-source streamfunction influence integral over a panel, evaluated at
# a field point. Includes the mfoil branch-cut shift.
function mfoil_constant_source_psi(xstar, ystar, logr1, logr2, lj, β)
    θ1 = atan(ystar, xstar)
    θ2 = atan(ystar, xstar - lj)
    if logr1 == 0
        θ1 = π
        θ2 = π
        β = zero(β)
    elseif logr2 == 0
        θ1 = zero(θ1)
        θ2 = zero(θ2)
        β = zero(β)
    end
    I3 = (ystar * (logr1 - logr2) - xstar * β + lj * θ2) / (2π)
    
    return (θ1 + θ2) > π ? (I3 - 0.25 * lj) : (I3 + 0.75 * lj)
end

# Linear-source streamfunction influence integrals over a panel. Returns
# (a, b) coefficients such that contribution is `a·σ_left + b·σ_right`.
function mfoil_linear_source_psi(xstar, ystar, logr1, logr2, lj, β, r1, r2)
    θ1 = atan(ystar, xstar)
    θ2 = atan(ystar, xstar - lj)

    # Branch-cut shifts for the linear-source integrals (mfoil's `linear=true`)
    if logr1 == 0
        β_eff = zero(β)
        θ2_eff = π
        θ1_eff = π
    elseif logr2 == 0
        β_eff = zero(β)
        θ2_eff = zero(θ2)
        θ1_eff = zero(θ1)
    else
        β_eff = β
        θ2_eff = θ2 < 0 ? θ2 + 2π : θ2
        θ1_eff = θ1 < 0 ? θ1 + 2π : θ1
    end

    I3 = (ystar * (logr1 - logr2) - xstar * β_eff + lj * θ2_eff) / (2π)
    I5 = (r2^2 * θ2_eff - r1^2 * θ1_eff - ystar * lj) / (4π)

    a = I3 - (xstar * I3 + I5) / lj
    b = (xstar * I3 + I5) / lj
    
    return a, b
end

# Linear-vortex velocity influence at a field point from a panel. Returns
# `(a, b)` where each is a `(Vx, Vy)` 2-tuple. Contribution to velocity is
# `a · γ_left + b · γ_right`.
function mfoil_linear_vortex_velocity(xstar, ystar, logr1, logr2, lj, β, panel)
    tx = panel.ctheta
    ty = panel.stheta
    nx = -ty
    ny = tx

    dlogr = logr1 - logr2
    Q1 = β / (2π)
    Q2 = (ystar * dlogr - xstar * β) / (2π * lj)
    u1 = Q1 + Q2
    u2 = -Q2

    Q3 = -dlogr / (2π)
    Q4 = (xstar * dlogr - lj + ystar * β) / (2π * lj)
    v1 = Q3 + Q4
    v2 = -Q4

    return (u1 * tx + v1 * nx, u1 * ty + v1 * ny),
           (u2 * tx + v2 * nx, u2 * ty + v2 * ny)
end

# Constant-source velocity influence at a field point from a panel.
function mfoil_constant_source_velocity(xstar, ystar, logr1, logr2, lj, panel)
    tx = panel.ctheta
    ty = panel.stheta
    nx = -ty
    ny = tx

    θ1 = atan(ystar, xstar)
    θ2 = atan(ystar, xstar - lj)
    u = (logr1 - logr2) / (2π)
    v = (θ2 - θ1) / (2π)

    return u * tx + v * nx, u * ty + v * ny
end

# Linear-source velocity influence at a field point from a panel. Returns
# `(a, b)` where each is a `(Vx, Vy)` 2-tuple. Contribution to velocity is
# `a · σ_left + b · σ_right`.
function mfoil_linear_source_velocity(xstar, ystar, logr1, logr2, lj, β, panel)
    tx = panel.ctheta
    ty = panel.stheta
    nx = -ty
    ny = tx

    dlogr = logr1 - logr2
    Q1 = dlogr / (2π)
    Q2 = (xstar * dlogr - lj + ystar * β) / (2π * lj)
    u1 = Q1 - Q2
    u2 = Q2

    Q3 = β / (2π)
    Q4 = (xstar * β - ystar * dlogr) / (2π * lj)
    v1 = Q3 - Q4
    v2 = Q4

    return (u1 * tx + v1 * nx, u1 * ty + v1 * ny),
           (u2 * tx + v2 * nx, u2 * ty + v2 * ny)
end
