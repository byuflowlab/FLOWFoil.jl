#=
Local Newton solver, single-station marchers, and the full boundary-layer
initialization march.

The local Newton (`local_newton`) is a 4-state Newton with two modes:
- **Direct mode** (initial): downstream `(θ, δ*, state3)` are the unknowns,
  ue is held fixed at F.ues[2]. The Jacobian is 3×3.
- **Inverse mode** (fallback): all four `(θ, δ*, state3, ue)` are unknowns,
  with an Hk-prescription equation augmenting the system to 4×4. Used
  when the direct mode encounters separation (Hk > Hkmax) or hits the
  iteration cap.

`station_laminar`, `station_turbulent`, `station_stagnation`, and
`station_transition` are thin wrappers that build a station NamedTuple
and dispatch to `local_newton`.

`initialize_boundary_layer` marches through the BL nodes once, starting
from a stagnation similarity solution, transitioning when either the
amplification factor crosses `ncrit` (free transition) or the arc-length
crosses `sft` (forced transition).
=#

# ── Local Newton ─────────────────────────────────────────────────────────────

"""
    local_newton(station, x0; tol=1e-10, max_iter=30, under_relaxation=true, iters_direct=12)

Solve a single-station BL residual via Newton iteration. Returns
`(x, converged::Bool)` where `x` has the same length as `x0` (typically 4:
`θ, δ*, state3, ue`). Internally switches between a 3×3 direct Newton
(treating ue as fixed) and a 4×4 inverse Newton (prescribing Hk instead).
"""
function local_newton(
    station,
    x0;
    tol=1e-10,
    max_iter=30,
    under_relaxation=true,
    iters_direct=12,
)
    x = copy(x0)
    direct = true
    Hk_step = zero(eltype(x))
    ds = station.xs[2] - station.xs[1]

    # Wrappers for the residual and for Hk at the downstream node.
    F(states) = station_residual(states, station)
    Hk_wrap(states) =
        get_Hk(states[2], states[1], station.ues[2], station.thermo;
               airfoil=station.airfoil, turbulent=station.turbulent,
               limit_Hk_max=false, limit_Hk_min=false)

    for i in 1:max_iter
        Fx = F(x)
        if norm(Fx) < tol
            return x, true
        end

        if direct
            # 3×3 direct solve: treat ue (= x[4]) as fixed
            xtemp = x[1:3]
            J = ForwardDiff.jacobian(F, xtemp)
            Δx = pinv(J) * (-Fx)
        else
            # 4×4 inverse solve: residuals + Hk-prescription equation
            J1 = ForwardDiff.jacobian(F, x)
            J2 = ForwardDiff.gradient(Hk_wrap, x)
            Jnew = vcat(J1, J2')
            Hk = Hk_wrap(x)
            b = vcat(-Fx, Hk_step - Hk)
            Δx = pinv(Jnew) * b
        end

        # Under-relaxation: limit fractional change in any state to 0.3
        if under_relaxation
            dm = max(abs(Δx[2] / station.δstar0), abs(Δx[1] / station.θ0))
            if !direct
                dm = max(dm, abs(Δx[4] / station.ues[1]))
            end
            if station.turbulent
                dm = max(dm, abs(Δx[3] / station.state30))
            elseif direct
                dm = max(dm, abs(Δx[3] / 10))
            end
            omega = dm > 0.3 ? 0.3 / dm : 1.0
            Δx = Δx .* omega
        end

        # Apply the update; in direct mode, ue (x[4]) stays untouched.
        xnew = if station.similar || !direct
            x + Δx
        else
            vcat(x[1:3] + Δx, x[4])
        end

        # Clamp state3 (n or cτ) to physically-sensible bounds
        if station.turbulent || station.transition
            xnew[3] = clamp(xnew[3], 1e-7, 0.3)
        else
            xnew[3] = max(xnew[3], zero(xnew[3]))
        end

        # Separation check — switch to inverse mode if Hk hits its cap
        Hk = Hk_wrap(xnew)
        Hkmax = station.turbulent ? 2.5 : 3.8

        if direct && ((Hk > Hkmax) || i > iters_direct || any(xnew[1:2] .< 0))
            # Switch to inverse mode: prescribe Hk for the next iterate based on
            # the upstream Hk and a stretch factor along the panel.
            direct = false
            Hk_up = get_Hk(
                station.δstar0,
                station.θ0,
                station.ues[1],
                station.thermo;
                airfoil=station.airfoil,
                turbulent=station.turbulent,
                limit_Hk_max=false,
                limit_Hk_min=false,
            )
            Hkr = ds / station.θ0

            if station.airfoil
                Hk_step = station.turbulent ? Hk_up - 0.15 * Hkr : Hk_up + 0.03 * Hkr
                Hk_step = max(Hk_step, Hkmax)
            else
                # In the wake, integrate Hk's ODE forward (6 Newton steps)
                H2 = Hk_up
                for _ in 1:6
                    H2 = H2 - (H2 + 0.03 * Hkr * (H2 - 1)^3 - Hk_up) /
                              (1 + 0.09 * Hkr * (H2 - 1)^2)
                end
                Hk_step = max(H2, 1.01)
            end
        else
            x = xnew
        end
    end

    Fx = F(x)
    
    return x, norm(Fx) < tol
end

# ── Single-station marchers ──────────────────────────────────────────────────

"""
    station_laminar(deltastar0, theta0, n0, ues, xs, thermo, state_refs; upper=true)

Solve one laminar BL station. Returns `(δ*, θ, n, ue)` at the downstream
node. If the Newton fails to converge, extrapolate values from the
upstream node.
"""
function station_laminar(
    deltastar0, theta0, n0, ues, xs, thermo, state_refs; upper=true, max_iter=30, tol=1e-10
)
    initial_values = [theta0, deltastar0, n0, ues[2]]
    station = bl_station(;
        xs=xs,
        δstar0=deltastar0,
        θ0=theta0,
        state30=n0,
        ues=ues,
        wake_gap=zeros(2),
        thermo=thermo,
        state_refs=state_refs,
        upper=upper,
        turbulent=false,
        transition=false,
        similar=false,
        airfoil=true,
        forced_transition_flag=false,
    )
    solution, converged = local_newton(station, initial_values; tol=tol, max_iter=max_iter)
    if !converged
        solution = _airfoil_station_extrapolation(deltastar0, theta0, n0, ues[2], xs)
    end
    
    return solution[2], solution[1], solution[3], solution[4]
end

"""
    station_turbulent(deltastar0, theta0, ct120, ues, xs, thermo, state_refs, wake_gap=zeros(2);
                       upper=true, airfoil=true)

Solve one turbulent BL station (airfoil or wake). Returns `(δ*, θ, cτ, ue)`.
"""
function station_turbulent(
    deltastar0,
    theta0,
    ct120,
    ues,
    xs,
    thermo,
    state_refs,
    wake_gap=zeros(2);
    upper=true,
    airfoil=true,
    max_iter=30,
    tol=1e-10,
)
    initial_values = [theta0, deltastar0, ct120, ues[2]]
    station = bl_station(;
        xs=xs,
        δstar0=deltastar0,
        θ0=theta0,
        state30=ct120,
        ues=ues,
        wake_gap=wake_gap,
        thermo=thermo,
        state_refs=state_refs,
        upper=upper,
        turbulent=true,
        transition=false,
        similar=false,
        airfoil=airfoil,
        forced_transition_flag=false,
    )
    solution, converged = local_newton(station, initial_values; tol=tol, max_iter=max_iter)
    if !converged
        solution = if airfoil
            _airfoil_station_extrapolation(deltastar0, theta0, ct120, ues[2], xs)
        else
            _wake_station_extrapolation(deltastar0, theta0, ct120, ues[2], xs)
        end
    end
    
    return solution[2], solution[1], solution[3], solution[4]
end

"""
    station_stagnation(deltastar0, theta0, n0, ues, xs, thermo, state_refs; upper=true)

Solve the stagnation similarity station. Returns the 3-state `(θ, δ*, n)`
at the stagnation point (ue is small but fixed and not solved for).
"""
function station_stagnation(
    deltastar0, theta0, n0, ues, xs, thermo, state_refs; upper=true, max_iter=30, tol=1e-10
)
    initial_values = [theta0, deltastar0, n0]
    station = bl_station(;
        xs=xs,
        δstar0=deltastar0,
        θ0=theta0,
        state30=n0,
        ues=ues,
        wake_gap=zeros(2),
        thermo=thermo,
        state_refs=state_refs,
        upper=upper,
        turbulent=false,
        transition=false,
        similar=true,
        airfoil=true,
        forced_transition_flag=false,
    )
    solution, converged = local_newton(station, initial_values; tol=tol, max_iter=max_iter)
    if !converged
        solution = _airfoil_station_extrapolation(deltastar0, theta0, n0, ues[2], xs)[1:3]
    end
    
    return solution
end

"""
    station_transition(deltastar1, theta1, n1, ues, xs, thermo, state_refs;
                        upper=true, forced_transition_flag=false)

Solve the laminar→turbulent transition station. Returns
`(δ*, θ, cτ, ue, xt, theta_t, delta_t)` at the downstream node and a few
auxiliaries (the last three are zero-valued placeholders; the actual transition 
location is stored in `state_refs.transition_upper / .transition_lower`).
"""
function station_transition(
    deltastar1,
    theta1,
    n1,
    ues,
    xs,
    thermo,
    state_refs;
    upper=true,
    forced_transition_flag=false,
    max_iter=30,
    tol=1e-10,
)
    ct2_0 = ct_at_transition(deltastar1, theta1, n1, ues[2], thermo; turbulent=false)
    initial_values = [theta1, deltastar1, ct2_0, ues[2]]
    station = bl_station(;
        xs=xs,
        δstar0=deltastar1,
        θ0=theta1,
        state30=n1,
        ues=ues,
        wake_gap=zeros(2),
        thermo=thermo,
        state_refs=state_refs,
        upper=upper,
        turbulent=false,
        transition=true,
        similar=false,
        airfoil=true,
        forced_transition_flag=forced_transition_flag,
    )
    solution, converged = local_newton(station, initial_values; tol=tol, max_iter=max_iter)
    if !converged
        solution = _airfoil_station_extrapolation(deltastar1, theta1, ct2_0, ues[2], xs)
    end
    
    return solution[2], solution[1], solution[3], solution[4], 0.0, 0.0, 0.0
end

# Extrapolation fallbacks when the local Newton doesn't converge.
function _airfoil_station_extrapolation(δstar0, θ0, state30, ue0, xs)
    return [θ0 * sqrt(xs[2] / xs[1]), δstar0 * sqrt(xs[2] / xs[1]), state30, ue0]
end

function _wake_station_extrapolation(δstar0, θ0, state30, ue0, xs)
    rlen = (xs[2] - xs[1]) / (10 * θ0)
    return [(θ0 + rlen * δstar0) / (1 + rlen), δstar0, state30, ue0]
end

# ── Stagnation interpolation ─────────────────────────────────────────────────

"""
    interpolate_stagnation_state(U1, U2, x1, x2)

Linear extrapolation of the BL state to the stagnation point at a small
positive x (`xst = 1e-6`), using the values at the first two BL nodes.
Edge velocity uses a quadratic extrapolation (ue → 0 at the stagnation
point linearly in arclength). Returns `(Ust, xst)`.
"""
function interpolate_stagnation_state(U1, U2, x1, x2)
    dx = x2 - x1
    rx = x2 / x1
    w1 = x2 / dx
    w2 = -x1 / dx

    Ust = U1 .* w1 .+ U2 .* w2

    # Quadratic-in-x extrapolation of ue
    wk1 = rx / dx
    wk2 = -1 / (rx * dx)
    K = wk1 * U1[4] + wk2 * U2[4]
    xst = eltype(U1)(1e-6)
    Ust[4] = K * xst

    # Sanity-clamp the state3 (n or cτ) to ≥ 0
    if Ust[3] < 0
        if Ust[3] isa ForwardDiff.Dual
            Ust[3] -= Ust[3].value
        else
            Ust[3] = 0
        end
    end

    return Ust, xst
end

# ── Wake first-node initialization ───────────────────────────────────────────

"""
    initialize_wake_first_node(UTE_l, UTE_u, Uw, te_geometry,
                                te_l_turbulent, te_u_turbulent, thermo, state_refs;
                                init_flag=false)

Bridge the boundary layer from the airfoil TE into the wake. The wake
starts at a single point just behind the TE midpoint; its initial state
is the sum of the upper- and lower-surface TE BL contributions plus the
TE gap thickness for δ*. If `init_flag=true`, returns the wake state
vector. Otherwise returns the residual `Uw - U_wake_first` for use as
part of the global Newton coupling.
"""
function initialize_wake_first_node(
    UTE_l,
    UTE_u,
    Uw,
    te_geometry,
    te_l_turbulent,
    te_u_turbulent,
    thermo,
    state_refs;
    init_flag=false,
)
    turbulent = true
    ctau_l = te_l_turbulent ? UTE_l[3] : ct_at_transition(UTE_l, thermo; turbulent)
    ctau_u = te_u_turbulent ? UTE_u[3] : ct_at_transition(UTE_u, thermo; turbulent)

    theta_l = UTE_l[1]
    theta_u = UTE_u[1]
    deltastar_l = UTE_l[2]
    deltastar_u = UTE_u[2]

    deltastar_wake = deltastar_l + deltastar_u + te_geometry.hTE
    theta_wake = theta_l + theta_u
    ctau_wake = (theta_l * ctau_l + theta_u * ctau_u) / (theta_l + theta_u)

    if init_flag
        return [theta_wake, deltastar_wake, ctau_wake]
    else
        return [Uw[1] - theta_wake, Uw[2] - deltastar_wake, Uw[3] - ctau_wake]
    end
end

# ── Full BL march ────────────────────────────────────────────────────────────

"""
    initialize_boundary_layer(s, thermo, state_refs, ues=nothing, U_stag_wake=nothing;
                               upper=true, airfoil=true, sft=nothing, max_iter=30, tol=1e-10)

March the boundary layer through the nodes of a single surface (airfoil
upper, airfoil lower, or wake). Starts at the stagnation point (airfoil
case) or at the wake-start state (wake case), then steps to each
downstream node by calling `station_*` at each step.

Free transition (laminar → turbulent) is triggered when the amplification
state crosses `ncrit`. Forced transition is triggered when `s[i+1] >
sft`. On detection, a transition station is solved at the same step, the
flow flag flips to turbulent, and the march continues.

# Returns
NamedTuple with `deltas, thetas, state3s, ues, idx_transition, st,
theta_t, delta_t`. `idx_transition` is the node index at which the flow
becomes turbulent (or `npanels + 2` if no transition).
"""
function initialize_boundary_layer(
    s,
    thermo,
    state_refs,
    ues=nothing,
    U_stag_wake=nothing;
    upper=true,
    airfoil=true,
    sft=nothing,
    max_iter=30,
    tol=1e-10,
)
    npanels = length(s) - 1
    idx_transition = npanels + 2
    forced_transition_flag = false

    deltas = zeros(npanels + 1)
    thetas = zeros(npanels + 1)
    state3s = zeros(npanels + 1)
    if isnothing(ues)
        ues = ones(npanels + 1)
    end

    if airfoil && isnothing(sft)
        sft = upper ? state_refs.sft_u[] : state_refs.sft_l[]
    end

    transition = false
    theta_t = 0.0
    delta_t = 0.0
    st = airfoil && forced_transition_flag ? sft : s[end]

    if airfoil
        K = s[1] > 1e-8 ? ues[1] / s[1] : ues[2] / s[2]
        θ0 = sqrt(0.45 * thermo.mu0 / (6 * thermo.rho0 * K))
        δstar0 = 2.2 * θ0
        ñ0 = 0.0

        turbulent = false
        hitstag = s[1] < 1e-8 * s[end]

        sstag = 1e-6
        ue_stag = hitstag ? ues[2] / s[2] * sstag : ues[1] / s[1] * sstag

        U_stag = station_stagnation(
            δstar0, θ0, ñ0, [ue_stag, ue_stag], [sstag, sstag], thermo, state_refs;
            upper=upper, max_iter=max_iter, tol=tol,
        )

        if hitstag
            thetas[1:2] .= U_stag[1]
            deltas[1:2] .= U_stag[2]
            state3s[1:2] .= U_stag[3]
            i = 2
        else
            thetas[1] = U_stag[1]
            deltas[1] = U_stag[2]
            state3s[1] = U_stag[3]
            i = 1
        end

        if airfoil && sft < s[2]
            turbulent = true
        end

        wake_gap = zeros(npanels + 1)
    else
        @assert !isnothing(U_stag_wake) "Wake initialization requires U_stag_wake."
        turbulent = true
        thetas[1] = U_stag_wake[1]
        deltas[1] = U_stag_wake[2]
        state3s[1] = U_stag_wake[3]
        i = 1
        wake_gap = state_refs.wake_gap
        sft = s[end]
    end

    while i <= npanels
        if airfoil && !turbulent && s[i + 1] > sft
            transition = true
            turbulent = true
            forced_transition_flag = true
            idx_transition = i + 1
        end

        if !turbulent
            δstar_temp, θ_temp, state3_temp, ue_temp = station_laminar(
                deltas[i],
                thetas[i],
                state3s[i],
                ues[i:(i + 1)],
                s[i:(i + 1)],
                thermo,
                state_refs;
                upper=upper, max_iter=max_iter, tol=tol,
            )
        elseif transition
            δstar_temp, θ_temp, state3_temp, ue_temp, st_new, theta_t, delta_t = station_transition(
                deltas[i],
                thetas[i],
                state3s[i],
                ues[i:(i + 1)],
                s[i:(i + 1)],
                thermo,
                state_refs;
                upper=upper,
                forced_transition_flag=forced_transition_flag,
                max_iter=max_iter, tol=tol,
            )
            st = st_new
            transition = false
            forced_transition_flag = false
        else
            δstar_temp, θ_temp, state3_temp, ue_temp = station_turbulent(
                deltas[i],
                thetas[i],
                state3s[i],
                ues[i:(i + 1)],
                s[i:(i + 1)],
                thermo,
                state_refs,
                wake_gap[i:(i + 1)];
                upper=upper, airfoil=airfoil, max_iter=max_iter, tol=tol,
            )
        end

        # Free-transition check: n exceeded ncrit during this laminar step.
        if !turbulent && state3_temp > thermo.ncrit
            transition = true
            turbulent = true
            idx_transition = i + 1
            i -= 1   # repeat this node as a transition station next iteration
        else
            deltas[i + 1] = δstar_temp
            thetas[i + 1] = θ_temp
            state3s[i + 1] = state3_temp
            ues[i + 1] = ue_temp
        end

        i += 1
    end

    return (;
        deltas=deltas,
        thetas=thetas,
        state3s=state3s,
        ues=ues,
        idx_transition=idx_transition,
        st=st,
        theta_t=theta_t,
        delta_t=delta_t,
    )
end
