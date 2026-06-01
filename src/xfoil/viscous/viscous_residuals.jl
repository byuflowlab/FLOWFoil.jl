#=
Boundary-layer station residual equations and the per-station residual
"closure" object used by the local Newton solver.

A *station* is a pair of streamwise-adjacent BL nodes; the residuals at
the downstream node (the one being solved for) depend on the upstream
state (fixed during the local Newton) and the downstream state (the
unknowns). Three equations per node:

1. Integrated momentum equation
2. Kinetic-energy / shape parameter equation
3. Amplification (laminar) or shear-lag (turbulent) equation

For a transition station, the integration is split into a laminar segment
upstream of the transition point and a turbulent segment downstream; the
two segments contribute additively to a single 3-component residual.

`bl_station(...)` builds a NamedTuple holding the fixed data; the
downstream-state residual is then computed by
`station_residual(states, station)`.
=#

"""
    bl_station(; xs, δstar0, θ0, state30, ues, wake_gap, thermo,
                  state_refs, upper, turbulent, transition, similar,
                  airfoil, forced_transition_flag)

Bundle all fixed data needed to evaluate the BL station residual.

# Arguments
- `xs`            : 2-element arc-length coordinates at the upstream / downstream nodes.
- `δstar0, θ0, state30` : upstream-node values held fixed during the local Newton.
- `ues`           : 2-element edge velocities at the two nodes.
- `wake_gap`      : 2-element wake-gap values at the two nodes (zero on the airfoil).
- `thermo`        : NamedTuple from `init_viscous_thermo`.
- `state_refs`    : NamedTuple of mutable scalar refs that the residual code may update:
                      `sft_l`, `sft_u`, `transition_lower`, `transition_upper`.
                    Each is a `Ref` so the global solver can read updated values back.
- `upper`         : true on upper surface, false on lower / wake.
- `turbulent`     : true if this station is in turbulent flow.
- `transition`    : true if this is a transition station (laminar→turbulent here).
- `similar`       : true if this is the stagnation station (similarity boundary condition).
- `airfoil`       : true on the airfoil, false in the wake.
- `forced_transition_flag` : true if transition has been forced at a specific x/c.
"""
function bl_station(;
    xs,
    δstar0,
    θ0,
    state30,
    ues,
    wake_gap=zeros(2),
    thermo,
    state_refs,
    upper=true,
    turbulent=false,
    transition=false,
    similar=false,
    airfoil=true,
    forced_transition_flag=false,
)
    return (;
        xs=xs,
        δstar0=δstar0,
        θ0=θ0,
        state30=state30,
        ues=ues,
        wake_gap=wake_gap,
        thermo=thermo,
        state_refs=state_refs,
        upper=upper,
        turbulent=turbulent,
        transition=transition,
        similar=similar,
        airfoil=airfoil,
        forced_transition_flag=forced_transition_flag,
    )
end

"""
    station_residual(states, station)

Evaluate the 3-component BL residual at the downstream node of a station.

`states` is one of:
- 3-vector: `(θ, δ*, state3)` for a station with prescribed edge velocity.
- 4-vector: `(θ, δ*, state3, ue)` for the inverse mode where ue is also an
  unknown (used in the local Newton when the direct mode hits a separation).

`station3` is the NamedTuple built by `bl_station`.
"""
function station_residual(states, station)
    if length(states) == 3
        if station.similar
            deltastars = [states[2], states[2]]
            thetas = [states[1], states[1]]
            state3s = [states[3], states[3]]
        else
            deltastars = [station.δstar0, states[2]]
            thetas = [station.θ0, states[1]]
            state3s = [station.state30, states[3]]
        end
        ues = station.ues
    elseif length(states) == 4
        if station.similar
            deltastars = [states[2], states[2]]
            thetas = [states[1], states[1]]
            state3s = [states[3], states[3]]
            ues = [states[4], states[4]]
        else
            deltastars = [station.δstar0, states[2]]
            thetas = [station.θ0, states[1]]
            state3s = [station.state30, states[3]]
            ues = [station.ues[1], states[4]]
        end
    else
        error("station_residual: states must have length 3 or 4 (got $(length(states))).")
    end

    if station.transition
        return residuals_transition(
            station.xs,
            deltastars,
            thetas,
            state3s,
            ues,
            station.thermo,
            station.state_refs;
            upper=station.upper,
            forced_transition_flag=station.forced_transition_flag,
        )
    else
        return bl_station_residuals(
            station.xs,
            deltastars,
            thetas,
            state3s,
            ues,
            station.thermo,
            station.wake_gap;
            airfoil=station.airfoil,
            turbulent=station.turbulent,
            upper=station.upper,
            similar=station.similar,
        )
    end
end

# ── Residuals (laminar / turbulent at a regular station) ─────────────────────

"""
    bl_station_residuals(xs, deltastars, thetas, state3s, ues, thermo, wake_gap=zeros(2);
                              airfoil=true, turbulent=true, upper=true, similar=false)

Standard BL station residual: 3-element vector of (momentum, shape,
amplification|shear-lag) residuals.
"""
function bl_station_residuals(
    xs,
    deltastars,
    thetas,
    state3s,
    ues,
    thermo,
    wake_gap=zeros(2);
    airfoil=true,
    turbulent=true,
    upper=true,
    similar=false,
)
    if !airfoil
        deltastars = deltastars .- wake_gap
    end

    Hwake1 = wake_gap[1] / thetas[1]
    Hwake2 = wake_gap[2] / thetas[2]
    Hwake = (Hwake1 + Hwake2) / 2

    H_vec = H_shape.(deltastars, thetas)
    M2 = get_Mach2.(ues, Ref(thermo))
    Reθ_vec = get_Reθ.(thetas, ues, Ref(thermo))

    ηup = upwind_factor(H_vec[1], H_vec[2], M2[1]; airfoil, turbulent)
    Hk_vec = get_Hk.(H_vec, M2; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)
    Hstarstar_vec = Hstarstar.(Hk_vec, M2)

    x_mid = (xs[1] + xs[2]) / 2
    δ_mid = (deltastars[1] + deltastars[2]) / 2
    θ_mid = (thetas[1] + thetas[2]) / 2
    H_mid = δ_mid / θ_mid
    ue_mid = (ues[1] + ues[2]) / 2
    M2_mid = get_Mach2(ue_mid, thermo)
    Hk_mid = get_Hk(H_mid, M2_mid; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)
    Reθ_mid = get_Reθ(θ_mid, ue_mid, thermo)

    if turbulent
        cτ12s = state3s

        if airfoil
            cf_vec = cf_turbulent.(Hk_vec, Reθ_vec, M2)
            cf_mid = cf_turbulent(Hk_mid, Reθ_mid, M2_mid)
        else
            TF = eltype(state3s)
            cf_vec = zeros(TF, 2)
            cf_mid = zero(TF)
        end
        cfxt_vec = cf_vec .* xs ./ thetas
        cfxt_mid = cf_mid * x_mid / θ_mid

        Hstar_vec = Hstar_turbulent.(Hk_vec, Reθ_vec, M2; airfoil)
        Us_vec = get_Us.(H_vec, Hk_vec, Hstar_vec; airfoil)
        cτeq12_vec = get_cτeq.(H_vec, Hk_vec, Hstar_vec, Reθ_vec, Us_vec; airfoil)
        cdi_vec = cdi_turbulent.(
            Hk_vec, Reθ_vec, cf_vec, Us_vec, Hstar_vec, cτ12s .^ 2; airfoil
        )
        cdixt_vec = cdi_vec .* xs ./ thetas
        δ_vec = bl_thickness.(deltastars, thetas, Hk_vec)

        R_mom = residual_momentum(
            xs, ηup, H_vec, ues, thetas, cfxt_vec, cfxt_mid, Hwake, thermo; airfoil
        )
        R_shape = residual_shape(
            xs,
            ηup,
            H_vec,
            ues,
            thetas,
            cfxt_vec,
            Hstar_vec,
            Hstarstar_vec,
            cdixt_vec,
            Hwake,
            thermo;
            airfoil,
        )
        R_lag = residual_shearlag(
            xs,
            ηup,
            H_vec,
            ues,
            deltastars,
            δ_vec,
            cτ12s,
            Us_vec,
            cτeq12_vec,
            cf_vec,
            Hk_vec,
            Reθ_vec,
            thermo;
            airfoil,
        )

        return [R_mom, R_shape, R_lag]
    else
        ñ_vec = state3s

        cf_vec = cf_laminar.(Hk_vec, Reθ_vec)
        cfxt_vec = cf_vec .* xs ./ thetas
        cf_mid = cf_laminar(Hk_mid, Reθ_mid)
        cfxt_mid = cf_mid * x_mid / θ_mid

        Hstar_vec = Hstar_laminar.(Hk_vec)
        cdi_vec = cdi_laminar.(Hk_vec, Reθ_vec)
        cdixt_vec = cdi_vec .* xs ./ thetas

        R_mom = residual_momentum(
            xs, ηup, H_vec, ues, thetas, cfxt_vec, cfxt_mid, Hwake, thermo; airfoil, similar
        )
        R_shape = residual_shape(
            xs,
            ηup,
            H_vec,
            ues,
            thetas,
            cfxt_vec,
            Hstar_vec,
            Hstarstar_vec,
            cdixt_vec,
            Hwake,
            thermo;
            airfoil,
            similar,
        )
        R_amp = residual_amplification(
            xs, H_vec, ues, thetas, ñ_vec, Hk_vec, Reθ_vec, thermo.ncrit; similar
        )

        return [R_mom, R_shape, R_amp]
    end
end

# ── Individual residual components ───────────────────────────────────────────

function residual_momentum(
    xs, ηup, H_vec, ues, thetas, cfxt_vec, cfxt_mid, Hwake, thermo; airfoil=true, similar=false
)
    M2s = get_Mach2.(ues, Ref(thermo))
    M2 = thermo.Mach > 0 ? upwind_avg(0.5, M2s[1], M2s[2]) : zero(eltype(ues))
    H = upwind_avg(0.5, H_vec[1], H_vec[2])

    # mfoil's Simpson-like average for the cfx/θ term (better drag accuracy)
    cfxt_avg = 0.25 * cfxt_vec[1] + 0.25 * cfxt_vec[2] + 0.5 * cfxt_mid

    uks = get_uk.(ues, Ref(thermo))

    if similar
        logθ = zero(eltype(uks))
        logue = one(eltype(uks))
        logx = one(eltype(uks))
    else
        logθ = log(thetas[2] / thetas[1])
        logue = log(uks[2] / uks[1])
        logx = log(xs[2] / xs[1])
    end

    return logθ + logue * (H + 2 + Hwake - M2) - logx * cfxt_avg / 2
end

function residual_shape(
    xs,
    ηup,
    H_vec,
    ues,
    θ_vec,
    cfxt_vec,
    Hstar_vec,
    Hstarstar_vec,
    cdixt_vec,
    Hwake,
    thermo;
    airfoil=true,
    similar=false,
)
    uks = get_uk.(ues, Ref(thermo))

    H = upwind_avg(0.5, H_vec[1], H_vec[2])
    Hstar = upwind_avg(0.5, Hstar_vec[1], Hstar_vec[2])
    Hstarstar = upwind_avg(ηup, Hstarstar_vec[1], Hstarstar_vec[2])
    cfxt = upwind_avg(ηup, cfxt_vec[1], cfxt_vec[2])
    cdixt = upwind_avg(ηup, cdixt_vec[1], cdixt_vec[2])

    if similar
        logHstar = zero(eltype(uks))
        logue = one(eltype(uks))
        logx = one(eltype(uks))
    else
        logHstar = log(Hstar_vec[2] / Hstar_vec[1])
        logue = log(uks[2] / uks[1])
        logx = log(xs[2] / xs[1])
    end

    return logHstar + logue * (2 * Hstarstar / Hstar + 1 - H - Hwake) +
           logx * (cfxt / 2 - cdixt)
end

function residual_amplification(xs, H_vec, ues, thetas, ns, Hk_vec, Reθ_vec, ncrit; similar=false)
    if similar
        return ns[2] + ns[1]   # stagnation: state3 = 0 at both nodes
    end
    dndx_vec = amplification_rate.(Hk_vec, Reθ_vec, ns, thetas; ncrit=ncrit)
    dndx = upwind_avg(0.5, dndx_vec[1], dndx_vec[2])
    
    return ns[2] - ns[1] - dndx * (xs[2] - xs[1])
end

function residual_amplification(xs, U1, U2, thermo)
    dndx = upwind_avg(0.5, amplification_rate(U1, thermo), amplification_rate(U2, thermo))
    
    return U2[3] - U1[3] - dndx * (xs[2] - xs[1])
end

function residual_shearlag(
    xs,
    ηup,
    H_vec,
    ues,
    deltastars,
    δ_vec,
    cτ12_vec,
    Us_vec,
    cτeq12_vec,
    cf_vec,
    Hk_vec,
    Reθ_vec,
    thermo;
    Klag=5.6,
    GB=0.75,
    airfoil=true,
)
    uks = get_uk.(ues, Ref(thermo))
    ηD = airfoil ? 1.0 : 0.9
    Δx = xs[2] - xs[1]

    δstar_avg = (deltastars[1] + deltastars[2]) / 2
    δ = upwind_avg(0.5, δ_vec[1], δ_vec[2])
    Us = upwind_avg(0.5, Us_vec[1], Us_vec[2])
    cτ12 = upwind_avg(ηup, cτ12_vec[1], cτ12_vec[2])
    cτeq12 = upwind_avg(ηup, cτeq12_vec[1], cτeq12_vec[2])
    cf = upwind_avg(ηup, cf_vec[1], cf_vec[2])
    Hk = upwind_avg(ηup, Hk_vec[1], Hk_vec[2])
    Reθ = upwind_avg(0.5, Reθ_vec[1], Reθ_vec[2])

    uq_val = uq(Hk, cf, Reθ, δstar_avg; airfoil)

    R = 2 * δ * log(cτ12_vec[2] / cτ12_vec[1]) -
        Klag / (GB * (1 + Us)) * Δx * (cτeq12 - ηD * cτ12) -
        2 * δ * (uq_val * Δx - log(uks[2] / uks[1]))
    
    return -R
end

# ── Transition station ───────────────────────────────────────────────────────

"""
    residuals_transition(xs, deltas, thetas, state3s, ues, thermo, state_refs;
                          upper=true, forced_transition_flag=false)

Residual at a transition station. The transition x is solved for via an
implicit ImplicitAD wrap of Brent's method (free transition) or read from
`state_refs.sft_u/sft_l` (forced transition). The laminar segment from
xs[1]→xt and the turbulent segment from xt→xs[2] each produce a 3-vector
of residuals; the returned residual is their sum (with the amplification
component zeroed when transition is forced).
"""
function residuals_transition(
    xs, deltas, thetas, state3s, ues, thermo, state_refs;
    upper=true, forced_transition_flag=false
)
    U1 = [thetas[1], deltas[1], state3s[1], ues[1]]
    U2 = [thetas[2], deltas[2], state3s[2], ues[2]]
    x1 = xs[1]
    dx = xs[2] - xs[1]

    Us = vcat(U1, U2, x1, dx)
    p = (; thermo)

    if forced_transition_flag
        xt = upper ? state_refs.sft_u[] : state_refs.sft_l[]
    else
        xt = ImplicitAD.implicit(_solve_for_xt_implicit, _xt_residual_implicit, Us, p)
    end

    # Record the transition location (real-valued part only)
    xt_real = eltype(xt) == Float64 ? xt : xt.value
    if upper
        state_refs.transition_upper[] = xt_real
    else
        state_refs.transition_lower[] = xt_real
    end

    w_turb = (xt - x1) / dx
    w_lam = 1 - w_turb
    Ut = w_lam .* U1 .+ w_turb .* U2
    Utl = copy(Ut)
    Utt = copy(Ut)
    if !forced_transition_flag
        Utl[3] = thermo.ncrit
    end
    Utt[3] = ct_at_transition(Utl, thermo; turbulent=true)

    x_lam = [xs[1], xt]
    x_turb = [xt, xs[2]]

    R_lam = bl_station_residuals(
        x_lam,
        [U1[2], Utl[2]],
        [U1[1], Utl[1]],
        [U1[3], Utl[3]],
        [U1[4], Utl[4]],
        thermo;
        airfoil=true,
        turbulent=false,
        upper=upper,
    )
    R_turb = bl_station_residuals(
        x_turb,
        [Utt[2], U2[2]],
        [Utt[1], U2[1]],
        [Utt[3], U2[3]],
        [Utt[4], U2[4]],
        thermo;
        airfoil=true,
        turbulent=true,
        upper=upper,
    )

    if forced_transition_flag
        R_lam[3] = 0
    end

    return R_lam .+ R_turb
end

# Solve for the transition x via Brent's method, with a Newton fallback.
function _solve_for_xt_implicit(Us, p)
    U1 = Us[1:4]
    U2 = Us[5:8]
    x1 = Us[9]
    dx = Us[10]
    p2 = (; p..., U1, U2, x1, dx)

    xt, info = FLOWMath.brent(s -> _xt_residual_scalar(s, p2), x1, x1 + dx; atol=1e-12)
    if info.flag != "CONVERGED"
        # Newton fallback from midpoint; if Newton also fails to converge,
        # default to the midpoint (mfoil's behavior).
        xt0 = x1 + 0.5 * dx
        xt_newton, converged = _newton_for_xt(xt0, p2)
        xt = converged ? xt_newton : xt0
    end
    
    return xt
end

function _xt_residual_implicit(y, Us, p)
    U1 = Us[1:4]
    U2 = Us[5:8]
    x1 = Us[9]
    dx = Us[10]
    p2 = (; p..., U1, U2, x1, dx)
    
    return _xt_residual_scalar(y[1], p2)
end

# Residual for the transition x: the integrated amplification from xs[1] to xt
# (using the trapezoidal rule for dn/dx) should equal ncrit - n_1.
function _xt_residual_scalar(xt, p)
    (; thermo, U1, U2, x1, dx) = p
    w_turb = (xt - x1) / dx
    w_lam = 1 - w_turb
    Ut = w_lam .* U1 .+ w_turb .* U2
    Ut = [Ut[1], Ut[2], thermo.ncrit, Ut[4]]
    dndx1 = amplification_rate(U1, thermo)
    dndxt = amplification_rate(Ut, thermo)
    
    return thermo.ncrit - U1[3] - (xt - x1) * (dndx1 + dndxt) / 2
end

# 1-D Newton on the transition x. Used as the Brent fallback. Returns
# `(xt, converged)`; the caller decides whether to use `xt` or default
# to the midpoint when `converged` is false. Step size tapers with `i`
# (mfoil convention) to encourage convergence near the bracket interior.
function _newton_for_xt(x0, p; tol=1e-12, max_iter=20)
    x = x0
    for i in 1:max_iter
        R = _xt_residual_scalar(x, p)
        if abs(R) < tol
            return x, true
        end
        dR = ForwardDiff.derivative(s -> _xt_residual_scalar(s, p), x)
        dx = -R / dR
        dmax = 0.2 * p.dx * (1.1 - i / max_iter)
        if abs(dx) > dmax
            dx = dx * dmax / abs(dx)
        end
        x += dx
    end
    R = _xt_residual_scalar(x, p)
    
    return x, abs(R) < tol
end

# ── Global residual: velocity equation ───────────────────────────────────────

"""
    velocity_residual(deltastars, ues, D, ueinv)

Velocity equation residual (`ue − (ueinv + D · (δ*·ue))`) used by the
global Newton coupling. All inputs are length `npanels_total + 1`.
"""
velocity_residual(deltastars, ues, D, ueinv) =
    ues .- (ueinv .+ D * (deltastars .* ues))
