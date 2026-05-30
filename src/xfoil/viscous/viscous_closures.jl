#=
Boundary-layer closure relations from Drela's xfoil/mfoil.

Where a closure needs thermodynamic state (Mach, KTb, KTl, H0, etc.) it
accepts a `thermo` NamedTuple built by `init_viscous_thermo`. Pure
incompressible-form closures take the relevant scalars directly.
=#

"""
    init_viscous_thermo(Re, Mach, rho, chord, Vinf; ncrit=9.0, GA=6.7, GB=0.75, GC=18.0, Cτ0=1.8, Eτ0=3.3)

Initialize the thermodynamic + closure-constants NamedTuple consumed by the
boundary-layer closures and residual functions. The closure constants
(ncrit, GA/GB/GC, Cτ0, Eτ0) are mfoil defaults.

For `Mach == 0` (the only fully-tested path), Karman-Tsien factors
degenerate to KTb=1, KTl=0 and stagnation quantities equal freestream
quantities.

# Returns
NamedTuple with: `Re, Mach, rho, mu, Vinf, chord, KTb, KTl, H0, mu0,
rho0, γair, ncrit, GA, GB, GC, Cτ0, Eτ0`.
"""
function init_viscous_thermo(
    Re, Mach, rho, chord, Vinf; ncrit=9.0, GA=6.7, GB=0.75, GC=18.0, Cτ0=1.8, Eτ0=3.3
)
    mu = rho * Vinf * chord / Re
    γair = 1.4
    γm1 = γair - 1
    KTb = 1.0
    KTl = 0.0
    H0 = 0.0
    mu0 = mu
    rho0 = rho
    Tsrat = 0.35
    if Mach > 0
        KTb = sqrt(1 - Mach^2)
        KTl = Mach^2 / (1 + KTb)^2
        H0 = (1 + 0.5 * γm1 * Mach^2) * Vinf^2 / (γm1 * Mach^2)
        Tr = 1 - 0.5 * Vinf^2 / H0
        finf = Tr^1.5 * (1 + Tsrat) / (Tr + Tsrat)
        mu0 = mu / finf
        rho0 = rho * (1 + 0.5 * γm1 * Mach^2)^(1 / γm1)
    end

    return (;
        Re,
        Mach,
        rho,
        mu,
        Vinf,
        chord,
        KTb,
        KTl,
        H0,
        mu0,
        rho0,
        γair,
        ncrit,
        GA,
        GB,
        GC,
        Cτ0,
        Eτ0,
    )
end

# ── Compressibility corrections ──────────────────────────────────────────────

"""
    get_uk(ue, thermo)

Karman-Tsien-corrected edge velocity. Identity when `Mach == 0`.
"""
function get_uk(ue, thermo)
    if thermo.Mach > 0
        den = 1 .- thermo.KTl .* (ue / thermo.Vinf) .^ 2
        return ue .* (1 - thermo.KTl) ./ den
    else
        return ue
    end
end

"""
    get_Mach2(ue, thermo)

Local Mach number squared.
"""
function get_Mach2(ue, thermo)
    if thermo.Mach > 0
        uk = get_uk(ue, thermo)
        γmi = thermo.γair - 1
        c2 = γmi * (thermo.H0 - 0.5 * uk^2)
        return uk^2 / c2
    else
        return zero(ue)
    end
end

"""
    get_density(ues, thermo)

Local get_density.
"""
function get_density(ues, thermo)
    if thermo.Mach > 0
        uk = get_uk(ues, thermo)
        M2 = get_Mach2(uk, thermo)
        γmi = thermo.γair - 1
        den = 1 + 0.5 * γmi * M2
        return thermo.rho0 / den^(1 / γmi)
    else
        return thermo.rho0
    end
end

"""
    get_Reθ(θ, ue, thermo)

Momentum-thickness Reynolds number.
"""
function get_Reθ(θ, ue, thermo)
    uk = get_uk(ue, thermo)
    if thermo.Mach > 0
        Ts = 0.35
        M2 = get_Mach2(ue, thermo)
        Tr = 1 - 0.5 * uk^2 / thermo.H0
        f = Tr^1.5 * (1 + Ts) / (Tr + Ts)
        mu1 = thermo.mu0 * f
        γmi = thermo.γair - 1
        den = 1 + 0.5 * γmi * M2
        rho1 = thermo.rho0 / den^(1 / γmi)
        return rho1 * uk * θ / mu1
    else
        return thermo.rho * uk * θ / thermo.mu
    end
end

# ── Shape parameters ─────────────────────────────────────────────────────────

"""
    H_shape(delta, theta)

Standard shape parameter `H = δ*/θ`. Specialized for `ForwardDiff.Dual`
to avoid divide-by-near-zero in derivative propagation.
"""
H_shape(delta, theta) = delta / theta

function H_shape(
    delta::ForwardDiff.Dual{Tag,T,N}, theta::ForwardDiff.Dual{Tag,T,N}
) where {Tag,T,N}
    H_val = delta.value / theta.value
    dH_ddelta = 1 / theta.value
    dH_dtheta = -H_val / theta.value
    partials = ntuple(
        i -> dH_ddelta * delta.partials[i] + dH_dtheta * theta.partials[i], N
    )
    
    
    return ForwardDiff.Dual{Tag}(H_val, partials)
end

function H_shape(delta::ForwardDiff.Dual{Tag,T,N}, theta::T) where {Tag,T,N}
    H_val = delta.value / theta
    dH_ddelta = 1 / theta
    partials = ntuple(i -> dH_ddelta * delta.partials[i], N)
    
    return ForwardDiff.Dual{Tag}(H_val, partials)
end

function H_shape(delta::T, theta::ForwardDiff.Dual{Tag,T,N}) where {Tag,T,N}
    H_val = delta / theta.value
    dH_dtheta = -H_val / theta.value
    partials = ntuple(i -> dH_dtheta * theta.partials[i], N)
    
    return ForwardDiff.Dual{Tag}(H_val, partials)
end

"""
    restrict_Hk(Hk; airfoil=true, turbulent=true, limit_Hk_min=false, limit_Hk_max=false)

Clamp the kinematic shape parameter `Hk` to numerically-sensible bounds.
"""
function restrict_Hk(
    Hk; airfoil=true, turbulent=true, limit_Hk_min=false, limit_Hk_max=false
)
    if limit_Hk_min
        Hk = airfoil ? max(Hk, 1.05) : max(Hk, 1.00005)
    end
    if limit_Hk_max
        Hk_cap = turbulent ? 2.5 : 3.8
        if Hk > Hk_cap
            Hk = Hk_cap
        end
    end
    
    return Hk
end

"""
    get_Hk(H, M2; airfoil=true, turbulent=true, limit_Hk_max=true, limit_Hk_min=true)

Kinematic shape parameter (compressibility-corrected).
"""
function get_Hk(
    H, M2; airfoil=true, turbulent=true, limit_Hk_max=true, limit_Hk_min=true
)
    Hk = (H - 0.29 * M2) / (1 + 0.113 * M2)
    
    return restrict_Hk(Hk; airfoil, turbulent, limit_Hk_max, limit_Hk_min)
end

get_Hk(δ, θ, thermo::NamedTuple; kwargs...) = get_Hk(δ / θ, thermo.Mach^2; kwargs...)

function get_Hk(δ, θ, ue, thermo::NamedTuple; kwargs...)
    return get_Hk(δ / θ, get_Mach2(ue, thermo); kwargs...)
end

"""
    bl_thickness(deltastar, theta, Hk)

Boundary-layer thickness.
"""
bl_thickness(deltastar, theta, Hk) = min(12 * theta, deltastar + theta * (3.15 + 1.72 / (Hk - 1)))

"""
    Hstarstar(Hk, M2)

Density shape parameter `H**`.
"""
Hstarstar(Hk, M2) = M2 * (0.064 / (Hk - 0.8) + 0.251)

# ── Kinetic-energy shape parameter ───────────────────────────────────────────

"""
    Hstar_laminar(Hk)

Kinetic-energy shape parameter `H* = δ_KE / θ` in laminar flow.
"""
function Hstar_laminar(Hk)
    Hk = restrict_Hk(Hk; airfoil=true, limit_Hk_min=true)
    H̃k = Hk - 4.35

    if Hk < 4.35
        return (0.0111 * H̃k^2 - 0.0278 * H̃k^3) / (Hk + 1) + 1.528 - 0.0002 * (H̃k * Hk)^2
    else
        return 0.015 * H̃k^2 / Hk + 1.528
    end
end

"""
    Hstar_turbulent(Hk, Reθ, M2; airfoil=true)

Kinetic-energy shape parameter in turbulent flow.
"""
function Hstar_turbulent(Hk, Reθ, M2; airfoil=true)
    Hk = restrict_Hk(Hk; airfoil, limit_Hk_min=true)
    Reθt = max(Reθ, 200)
    H0 = min(3 + 400 / Reθ, 4)
    AHstar = Hk - H0 + 4 / log(Reθt)
    Hr = (H0 - Hk) / (H0 - 1)
    if Hk < H0
        Hstar_incomp =
            1.5 + 4 / Reθt + 1.5 * Hr^2 / (Hk + 0.5) * (0.5 - 4 / Reθt)
    else
        Hstar_incomp =
            1.5 +
            4 / Reθt +
            (Hk - H0)^2 * (0.015 / Hk + 0.007 * log(Reθt) / AHstar^2)
    end
    
    return (Hstar_incomp + 0.028 * M2) / (1 + 0.014 * M2)
end

# Convenience overloads that take a state vector U = (θ, δ*, n_or_cτ, ue) and thermo.
function Hstar_laminar(δstar, θ, _n, ue, thermo)
    Hk = get_Hk(δstar / θ, get_Mach2(ue, thermo))
    
    return Hstar_laminar(Hk)
end
Hstar_laminar(U, thermo) = Hstar_laminar(U[2], U[1], U[3], U[4], thermo)

function Hstar_turbulent(δstar, θ, _ct, ue, thermo; kwargs...)
    M2 = get_Mach2(ue, thermo)
    Hk = get_Hk(δstar / θ, M2)
    Reθ = get_Reθ(θ, ue, thermo)
    
    return Hstar_turbulent(Hk, Reθ, M2; kwargs...)
end

# ── Friction coefficient ─────────────────────────────────────────────────────

"""
    cf_laminar(Hk, Reθ)

Skin-friction coefficient in laminar flow (divided by Reθ).
"""
function cf_laminar(Hk, Reθ)
    cf = if real(Hk) < 5.5
        0.0727 * (5.5 - Hk)^3 / (Hk + 1) - 0.07
    else
        0.015 * (1 - 1 / (Hk - 4.5))^2 - 0.07
    end

    return cf / Reθ
end

"""
    cf_wake()

Skin friction in the wake (always zero).
"""
cf_wake() = 0.0

"""
    cf_turbulent(Hk, Reθ, M2; γair=1.4)

Skin-friction coefficient in turbulent flow.
"""
function cf_turbulent(Hk, Reθ, M2; γair=1.4)
    Fc = sqrt(1 + 0.5 * (γair - 1) * M2)
    Acf = -1.33 * Hk
    if real(Acf) < -17
        Acf = -20 + 3 * exp((Acf + 17) / 3)
    end
    Bcf = max(log10(Reθ / Fc), 1.303)
    cf = 0.3 * exp(Acf) * Bcf^(-1.74 - 0.31 * Hk) +
         0.00011 * (tanh(4 - Hk / 0.875) - 1)
    
         return cf / Fc
end

# ── Dissipation coefficient ──────────────────────────────────────────────────
"""
    cdi_laminar(Hk, Reθ)

Dissipation coefficient in laminar flow.
"""
function cdi_laminar(Hk, Reθ)
    cdi = if real(Hk) < 4
        0.00205 * (4 - Hk)^5.5 + 0.207
    else
        -(0.0016 * (Hk - 4)^2) / (1 + 0.02 * (Hk - 4)^2) + 0.207
    end
    
    return cdi / Reθ
end

cdi_wall(Hk, Reθ, cf, Us, Hstar) =
    (cf * Us / (2 * Hstar)) * (1 + tanh(((Hk - 1) * log(Reθ)) / 2.1))

cdi_outer(cτ, Hstar, Us) = 2 * cτ * (0.995 - Us) / Hstar

cdi_stress(Reθ, Hstar, Us) = 0.3 * (0.995 - Us)^2 / (Hstar * Reθ)

function cdi_lamwake(Hk, Reθ)
    Hstar = Hstar_laminar(Hk)
    
    return 2.2 * (1 - 1 / Hk)^2 / (Hk * Hstar * Reθ)
end

"""
    cdi_laminar(Hk, Reθ)

Dissipation coefficient in laminar flow.
"""
function cdi_turbulent(Hk, Reθ, cf, Us, Hstar, cτ; airfoil=true)
    wall = airfoil ? cdi_wall(Hk, Reθ, cf, Us, Hstar) : 0.0
    outer = cdi_outer(cτ, Hstar, Us)
    stress = cdi_stress(Reθ, Hstar, Us)
    laminar = airfoil ? cdi_laminar(Hk, Reθ) : cdi_lamwake(Hk, Reθ)
    scale = airfoil ? 1.0 : 2.0
    
    return max(wall + outer + stress, laminar) * scale
end

cdi_wake(cdi_outer, cdi_stress, cdi_lamwake) =
    2 * max(cdi_outer + cdi_stress, cdi_lamwake)

# ── Shear-lag closures ───────────────────────────────────────────────────────

"""
    get_Us(H, Hk, Hstar; GB=0.75, airfoil=true)

Normalized wall slip velocity.
"""
function get_Us(H, Hk, Hstar; GB=0.75, airfoil=true)
    Hk = restrict_Hk(Hk; airfoil, limit_Hk_min=true)
    Us = Hstar / 2 * (1 - (Hk - 1) / (H * GB))
    
    return airfoil ? min(Us, 0.98) : min(Us, 0.99995)
end

function uq(Hk, cf, Reθ, δstar; GA=6.7, GB=0.75, airfoil=true)
    Hkc = get_Hkc(Hk, Reθ; airfoil)
    ηD = airfoil ? 1.0 : 0.9
    
    return (0.5 * cf - (Hkc / (GA * ηD * Hk))^2) / (GB * δstar)
end

function get_Hkc(Hk, Reθ; GC=18.0, airfoil=true)
    GC = airfoil ? GC : 0.0
    
    return max(Hk - 1 - GC / Reθ, 0.01)
end

"""
    get_cτeq(H, Hk, Hstar, Reθ, Us; GA=6.7, GB=0.75, airfoil=true)

Equilibrium value of sqrt(cτ) for the shear-lag equation.
"""
function get_cτeq(H, Hk, Hstar, Reθ, Us; GA=6.7, GB=0.75, airfoil=true)
    Hk = restrict_Hk(Hk; airfoil, limit_Hk_min=true)
    Hkc = get_Hkc(Hk, Reθ; airfoil)
    # Guard against tiny negative argument from floating-point roundoff
    # near separation; the value should be ≥ 0 mathematically.
    arg = Hstar * (Hk - 1) * Hkc^2 / (2 * GA^2 * GB * (1 - Us) * H * Hk^2)
    
    return sqrt(max(arg, zero(arg)))
end

"""
    ct_at_transition(delta, theta, n, ue, thermo; turbulent=false)

Initial sqrt(cτ) for a transition station (mfoil's transition shear-stress
initializer).
"""
function ct_at_transition(delta, theta, n, ue, thermo; turbulent=false)
    Cτ0 = thermo.Cτ0
    Eτ0 = thermo.Eτ0
    M2 = get_Mach2(ue, thermo)
    H = delta / theta
    Hk = get_Hk(H, M2; turbulent, limit_Hk_max=false)
    Reθ = get_Reθ(theta, ue, thermo)
    Hstar = turbulent ? Hstar_turbulent(Hk, Reθ, thermo.Mach; airfoil=true) :
                        Hstar_laminar(Hk)
    Us = get_Us(H, Hk, Hstar)
    cτeq = get_cτeq(H, Hk, Hstar, Reθ, Us)

    return Cτ0 * exp(-Eτ0 / (Hk - 1)) * cτeq
end

function ct_at_transition(U, thermo; turbulent=false)
    return ct_at_transition(U[2], U[1], U[3], U[4], thermo; turbulent)
end

# ── Laminar amplification ────────────────────────────────────────────────────

amplification_rate(rn, fn, gn, ϵn, θ) = (rn * fn * gn + ϵn) / θ

function amplification_rate(Hk, Reθ, n, θ; ncrit=9)
    Hk = restrict_Hk(Hk; airfoil=true, limit_Hk_min=true)
    Ĥ = get_Ĥ(Hk)
    L0 = get_L0(Ĥ)
    sn = get_sn(Reθ, L0)
    rn = get_rn(sn)
    fn = get_fn(Ĥ)
    gn = get_gn(Ĥ)
    ϵn = get_ϵn(n, ncrit)

    return amplification_rate(rn, fn, gn, ϵn, θ)
end

"""
    amplification_rate(U, thermo)

dndx based on states and thermo parameters.
"""
function amplification_rate(U, thermo)
    M2 = get_Mach2(U[4], thermo)
    H = U[2] / U[1]
    Hk = get_Hk(H, M2; turbulent=false, limit_Hk_min=true, limit_Hk_max=false)
    Reθ = get_Reθ(U[1], U[4], thermo)

    return amplification_rate(Hk, Reθ, U[3], U[1]; ncrit=thermo.ncrit)
end

get_fn(Ĥ) = -0.05 + 2.7 * Ĥ - 5.5 * Ĥ^2 + 3 * Ĥ^3 + 0.1 * exp(-20 * Ĥ)
get_gn(Ĥ) = 0.028 / Ĥ - 0.0345 * exp(-(3.87 * Ĥ - 2.52)^2)

function get_rn(sn)
    if real(sn) < 0
        return zero(sn)
    elseif real(sn) >= 1
        return one(sn)
    else
        return 3 * sn^2 - 2 * sn^3
    end
end

get_sn(Reθ, L0) = (log10(max(Reθ, eps(Reθ))) - (L0 - 0.1)) / 0.2
get_L0(Ĥ) = 2.492 * Ĥ^0.43 + 0.7 * (1 + tanh(14 * Ĥ - 9.24))
get_Ĥ(Hk) = 1 / (Hk - 1)
get_ϵn(n, ncrit) = 0.001 * (1 + tanh(5 * (n - ncrit)))

# ── Upwinding ────────────────────────────────────────────────────────────────

"""
    upwind_avg(ηup, x1, x2)

Linear interpolation factor used in midpoint-rule integration of the BL
equations. `ηup = 0.5` is centered; `ηup → 1` is fully upwind.
"""
upwind_avg(ηup, x1, x2) = (1 - ηup) * x1 + ηup * x2

"""
    upwind_factor(H1, H2, M2; airfoil=true, turbulent=true)

Upwinding factor used by the BL momentum/shape residuals. Smoothly biases
toward upwind when the downstream Hk diverges from the upstream Hk (mfoil
formula). `Cup = 5` on airfoils, `Cup = 1` in the wake.
"""
function upwind_factor(H1, H2, M2; airfoil=true, turbulent=true)
    Hk1 = get_Hk(H1, M2; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)
    Hk2 = get_Hk(H2, M2; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)
    fHu = (Hk2 - 1) / (Hk1 - 1)
    Cup = airfoil ? 5 : 1

    return 1 - 0.5 * exp(-log(abs(fHu))^2 * Cup / Hk2^2)
end
