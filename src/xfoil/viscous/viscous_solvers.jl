#=
Global coupled-viscous Newton solver and supporting helpers.

The coupled state vector has 4 components per node (`θ, δ*, state3, ue`)
on the airfoil + wake. The residual at each node is 4 entries: the three
station residuals (`R_mom, R_shape, R_amp/lag`) plus the velocity residual
(`ue = ueinv + D · (δ*·ue)`).

Per-iteration auxiliary updates:
- `stagnation_point_move!` shifts the stagnation-point node index when
  the Newton tries to flip the sign of `ue` on the first lower/upper node.
- `resolve_transition_location!` marches the amplification ODE again to
  re-detect the transition location given the current state.

`build_global_residuals` packages the per-problem fixed data (geometry,
influence matrices, mutable index vectors, etc.) into a NamedTuple; the
residual function `global_residual(U, gr)` consumes it.
=#

"""
    build_global_residuals(; s, su, sl, sw, idx_u, idx_l, idx_w,
                            idx_transition_u, idx_transition_l, turbulent_nodes,
                            thermo, state_refs, Bprime, Dwake, D, d, ueinv,
                            airfoil_panels, te_geometry, wake_panels, swake0,
                            Vinf, Mach, verbose)

Bundle the global-residual closure data. Several fields are *mutated* during
the Newton iterations: `idx_u`/`idx_l` (stagnation shifts), `D`/`d`/`ueinv`
(after stagnation moves), `idx_transition_u/_l`, `turbulent_nodes`,
`state_refs.transition_*`. Pass them as Vectors/Refs so the closure sees
the updates.
"""
function build_global_residuals(;
    s, su, sl, sw,
    idx_u, idx_l, idx_w,
    idx_transition_u, idx_transition_l,
    turbulent_nodes,
    thermo, state_refs,
    Bprime, Dwake, D, d, ueinv,
    airfoil_panels, te_geometry, wake_panels,
    swake0,
    Vinf, Mach,
    verbose=false,
)
    return (;
        s, su, sl, sw,
        idx_u, idx_l, idx_w,
        idx_transition_u, idx_transition_l,
        turbulent_nodes,
        thermo, state_refs,
        Bprime, Dwake, D, d, ueinv,
        airfoil_panels, te_geometry, wake_panels,
        swake0,
        Vinf, Mach,
        verbose,
    )
end

# ── Global residual evaluation ──────────────────────────────────────────────

"""
    global_residual(U, gr)

Residual vector at every node and the velocity equation residual at every
node. `U` is the flattened state vector `(θ_1, δ*_1, state3_1, ue_1,
θ_2, δ*_2, ...)` of length `4 · n_nodes_total`. Returns a flat vector of
length `4 · n_nodes_total` matching the layout.
"""
function global_residual(U, gr)
    U_nodes = reshape(U, 4, :)'
    nnodes = size(U_nodes, 1)
    Rs_global = zeros(eltype(U_nodes), nnodes, 3)

    # Re-split arc-length around the current stagnation point implied by the ue state
    su, sl, sw, _ = resplit_at_stagnation(gr.s, U_nodes[:, 4], gr.idx_u, gr.idx_l, gr.idx_w, gr.swake0)

    # ── Lower surface ───────────────────────────────────────────────────────
    if sl[1] < 1e-8 * sl[end]
        i0 = 2
        idx_stagnation = gr.idx_l[1] - 1
    else
        i0 = 1
        idx_stagnation = gr.idx_l[1]
    end

    Ust, xst = interpolate_stagnation_state(
        U_nodes[idx_stagnation, :], U_nodes[idx_stagnation - 1, :], sl[i0], sl[i0 + 1]
    )

    stag_station = bl_station(;
        xs=[xst, xst],
        δstar0=Ust[2], θ0=Ust[1], state30=Ust[3],
        ues=[Ust[4], Ust[4]],
        wake_gap=zeros(eltype(U), 2),
        thermo=gr.thermo, state_refs=gr.state_refs,
        upper=false, turbulent=false, transition=false, similar=true, airfoil=true,
        forced_transition_flag=false,
    )
    Rs_global[idx_stagnation, :] .= station_residual(Ust[1:3], stag_station)

    if i0 == 2
        Rs_global[gr.idx_l[1], :] = U_nodes[gr.idx_l[1], 1:3] .- Ust[1:3]
    end

    turbulent = false
    transition = false
    forced_transition_flag = false
    idx_surface_l = gr.idx_l[1]:-1:gr.idx_l[end]
    for i in (i0 + 1):length(idx_surface_l)
        idx = idx_surface_l[i]
        if !turbulent && (i == gr.idx_transition_l[1] || sl[i] > gr.state_refs.sft_l[])
            transition = true
            forced_transition_flag = sl[i] > gr.state_refs.sft_l[]
        end
        station = bl_station(;
            xs=eltype(U).(sl[(i - 1):i]),
            δstar0=U_nodes[idx + 1, 2], θ0=U_nodes[idx + 1, 1], state30=U_nodes[idx + 1, 3],
            ues=U_nodes[(idx + 1):-1:idx, 4],
            wake_gap=zeros(eltype(U), 2),
            thermo=gr.thermo, state_refs=gr.state_refs,
            upper=false, turbulent=turbulent, transition=transition, similar=false,
            airfoil=true, forced_transition_flag=forced_transition_flag,
        )
        Rs_global[idx, :] .= station_residual(U_nodes[idx, 1:3], station)
        if transition
            transition = false
            forced_transition_flag = false
            turbulent = true
        end
        if turbulent
            gr.turbulent_nodes[idx] = true
        end
    end
    TE_l_isturbulent = turbulent

    # ── Upper surface ───────────────────────────────────────────────────────
    if su[1] < 1e-8 * su[end]
        i0 = 2
        idx_stagnation = gr.idx_u[1] + 1
    else
        i0 = 1
        idx_stagnation = gr.idx_u[1]
    end

    Ust, xst = interpolate_stagnation_state(
        U_nodes[idx_stagnation, :], U_nodes[idx_stagnation + 1, :], su[i0], su[i0 + 1]
    )

    stag_station = bl_station(;
        xs=[xst, xst],
        δstar0=Ust[2], θ0=Ust[1], state30=Ust[3],
        ues=[Ust[4], Ust[4]],
        wake_gap=zeros(eltype(U), 2),
        thermo=gr.thermo, state_refs=gr.state_refs,
        upper=true, turbulent=false, transition=false, similar=true, airfoil=true,
        forced_transition_flag=false,
    )
    Rs_global[idx_stagnation, :] .= station_residual(Ust[1:3], stag_station)

    if i0 == 2
        Rs_global[gr.idx_u[1], :] = U_nodes[gr.idx_u[1], 1:3] .- Ust[1:3]
    end

    turbulent = false
    transition = false
    forced_transition_flag = false
    idx_surface_u = gr.idx_u[1]:gr.idx_u[end]
    for i in (i0 + 1):length(idx_surface_u)
        idx = idx_surface_u[i]
        if !turbulent && (i == gr.idx_transition_u[1] || su[i] > gr.state_refs.sft_u[])
            transition = true
            forced_transition_flag = su[i] > gr.state_refs.sft_u[]
        end
        station = bl_station(;
            xs=eltype(U).(su[(i - 1):i]),
            δstar0=U_nodes[idx - 1, 2], θ0=U_nodes[idx - 1, 1], state30=U_nodes[idx - 1, 3],
            ues=U_nodes[(idx - 1):idx, 4],
            wake_gap=zeros(eltype(U), 2),
            thermo=gr.thermo, state_refs=gr.state_refs,
            upper=true, turbulent=turbulent, transition=transition, similar=false,
            airfoil=true, forced_transition_flag=forced_transition_flag,
        )
        Rs_global[idx, :] .= station_residual(U_nodes[idx, 1:3], station)
        if transition
            transition = false
            forced_transition_flag = false
            turbulent = true
        end
        if turbulent
            gr.turbulent_nodes[idx] = true
        end
    end
    TE_u_isturbulent = turbulent

    # ── Wake ────────────────────────────────────────────────────────────────
    for i in 1:(length(gr.wake_panels) + 1)
        if i == 1
            Rs = initialize_wake_first_node(
                U_nodes[gr.idx_l[end], :],
                U_nodes[gr.idx_u[end], :],
                U_nodes[gr.idx_w[1], :],
                gr.te_geometry,
                TE_l_isturbulent,
                TE_u_isturbulent,
                gr.thermo,
                gr.state_refs;
                init_flag=false,
            )
            Rs_global[gr.idx_w[1], :] .= Rs
        else
            station = bl_station(;
                xs=eltype(U).(sw[(i - 1):i]),
                δstar0=U_nodes[gr.idx_w[i - 1], 2],
                θ0=U_nodes[gr.idx_w[i - 1], 1],
                state30=U_nodes[gr.idx_w[i - 1], 3],
                ues=U_nodes[gr.idx_w[i - 1]:gr.idx_w[i], 4],
                wake_gap=eltype(U).(gr.state_refs.wake_gap[(i - 1):i]),
                thermo=gr.thermo, state_refs=gr.state_refs,
                upper=true, turbulent=true, transition=false, similar=false,
                airfoil=false, forced_transition_flag=false,
            )
            Rs_global[gr.idx_w[i], :] .= station_residual(U_nodes[gr.idx_w[i], 1:3], station)
        end
    end

    # Velocity equation residual at every node
    Rs_velocity = velocity_residual(U_nodes[:, 2], U_nodes[:, 4], gr.D, gr.ueinv)

    return vcat(reshape(Rs_global', :, 1)[:, 1], Rs_velocity)
end

# ── Stagnation point bookkeeping ────────────────────────────────────────────

"""
    stagnation_point_move!(gr, ue)

If the Newton update has flipped the sign of `ue` at the first upper or
lower BL node, shift the stagnation point one or more panels to keep all
BL `ue` values positive. Flips the sign of `ue, d, ueinv` on the panels
spanned by the move and shifts `idx_u`, `idx_l`, `idx_transition_u/l` in
place. If a move happens, the `D` matrix is recomputed.

Returns `(su, sl, sw, s_stag, ue)`.
"""
function stagnation_point_move!(gr, ue)
    newpanel_flag = false

    if (ue[gr.idx_u[1]] > 0 && any(ue[gr.idx_u] .< 0)) ||
       (ue[gr.idx_l[1]] > 0 && any(ue[gr.idx_l] .< 0))
        error("Negative ue past the first node on a surface — Newton step is unphysical.")
    end

    if ue[gr.idx_u[1]] < 0
        # move stagnation up
        steps = findfirst(x -> x > 0, ue[gr.idx_u[1]:gr.idx_u[end]]) - 1
        newpanel_idx = gr.idx_u[1] + steps
        for i in gr.idx_u[1]:(newpanel_idx - 1)
            ue[i] = -ue[i]
            gr.d[i] = -gr.d[i]
            gr.ueinv[i] *= -1
        end
        gr.idx_u[1] += steps
        gr.idx_l[1] += steps
        gr.idx_transition_u[1] -= steps
        gr.idx_transition_l[1] += steps
        newpanel_flag = true
    elseif ue[gr.idx_l[1]] < 0
        # move stagnation down
        steps = findfirst(x -> x > 0, ue[gr.idx_l[1]:-1:gr.idx_l[end]]) - 1
        newpanel_idx = gr.idx_l[1] - steps
        for i in gr.idx_l[1]:-1:(newpanel_idx + 1)
            ue[i] = -ue[i]
            gr.d[i] = -gr.d[i]
            gr.ueinv[i] *= -1
        end
        gr.idx_u[1] -= steps
        gr.idx_l[1] -= steps
        gr.idx_transition_u[1] += steps
        gr.idx_transition_l[1] -= steps
        newpanel_flag = true
    end

    su, sl, sw, s_stag = resplit_at_stagnation(gr.s, ue, gr.idx_u, gr.idx_l, gr.idx_w, gr.swake0)

    if newpanel_flag
        # Recompute Dprime and D with the updated direction vector
        d_airfoil = gr.d[1:gr.idx_u[end]]
        Dprime_new = create_Dprime_influence_matrix(gr.s, gr.sw, d_airfoil, gr.d)
        Dnew = spdiagm(gr.d) * vcat(gr.Bprime, gr.Dwake) * Dprime_new
        gr.D[:] = Dnew[:]
    end

    return su, sl, sw, s_stag, ue
end

# ── Transition relocation (per-iteration) ───────────────────────────────────

"""
    resolve_transition_location!(gr, U; upper=true)

Re-march the amplification ODE along the surface to redetect the
transition node. If transition moved, fills in fresh turbulent `cτ`
states between the old and new transition nodes (or marks newly-laminar
nodes if transition moved downstream) and updates
`gr.idx_transition_u/_l` and `gr.turbulent_nodes` in place. `U` is the
`(nnodes_total × 4)` state matrix (mutated).
"""
function resolve_transition_location!(gr, U; upper=true)
    if upper
        idx_surface = gr.idx_u
        idx_nodes = idx_surface[1]:idx_surface[end]
        idx_transition = gr.idx_transition_u
    else
        idx_surface = gr.idx_l
        idx_nodes = idx_surface[1]:-1:idx_surface[end]
        idx_transition = gr.idx_transition_l
    end

    su, sl, _, _ = resplit_at_stagnation(gr.s, U[:, 4], gr.idx_u, gr.idx_l, gr.idx_w, gr.swake0)
    s = upper ? su : sl

    amp_states = U[idx_nodes, 3]   # save in case we revert

    idx_turb = march_amplification!(s, U, gr.turbulent_nodes, idx_nodes, gr; upper=upper)

    if idx_turb == idx_transition[1]
        # No change — restore amplification states (the march may have rewritten them)
        U[idx_nodes, 3] .= amp_states
        return
    end

    if gr.verbose
        surface = upper ? "upper" : "lower"
        println("Update $surface transition: last lam [$(idx_transition[1]-1)]->[$(idx_turb-1)]")
    end

    if idx_turb < idx_transition[1]
        # Transition moved upstream — fill in turbulent cτ between idx_turb and old transition
        ct0 = ct_at_transition(U[idx_nodes[idx_turb - 1], :], gr.thermo; turbulent=true)
        if idx_transition[1] < length(s)
            ct1 = U[idx_nodes[idx_transition[1]], 3]
        else
            ct1 = ct0
        end
        dx = s[min(idx_transition[1], length(s))] - s[idx_turb]

        for i in (idx_turb):min(idx_transition[1] - 1, length(idx_nodes))
            f = if dx == 0 || i == idx_turb
                0.0
            elseif idx_turb == idx_transition[1] - 1
                1.0
            else
                (s[i] - s[idx_turb]) / dx
            end
            U[idx_nodes[i], 3] = ct0 + f * (ct1 - ct0)
            gr.turbulent_nodes[idx_nodes[i]] = true
        end
    else
        # Transition moved downstream — unmark nodes between old and new transitions
        for i in idx_transition[1]:(idx_turb - 1)
            gr.turbulent_nodes[idx_nodes[i]] = false
        end
    end

    idx_transition[1] = idx_turb
    
    return
end

"""
    march_amplification!(s, U, turbulent_nodes, idx_nodes, gr; upper=true)

Re-march the laminar amplification equation along the surface to find
the first node where `n > ncrit` or `s > sft`. Returns the transition
node index (or `length(idx_nodes) + 1` if no transition).
"""
function march_amplification!(s, U, turbulent_nodes, idx_nodes, gr; upper=true)
    U[idx_nodes[1], 3] = 0
    sft = upper ? gr.state_refs.sft_u[] : gr.state_refs.sft_l[]

    for i in 2:length(idx_nodes)
        params_step = (;
            xs=[s[i - 1], s[i]],
            U1=U[idx_nodes[i - 1], :],
            U2=U[idx_nodes[i], :],
            thermo=gr.thermo,
            turbulent=false,
        )
        n0 = turbulent_nodes[idx_nodes[i]] ? U[idx_nodes[i - 1], 3] * 1.01 :
                                              U[idx_nodes[i], 3]
        solution, converged = local_newton_amp(_amplification_residual, [n0], params_step;
                                                tol=1e-12, max_iter=50)
        if !converged
            solution = [n0]
        end
        n2 = solution[1]

        if n2 > gr.thermo.ncrit || s[i] > sft
            return i
        end
        U[idx_nodes[i], 3] = n2
    end
    
    return length(idx_nodes) + 1
end

# 1-D residual for the amplification re-march: residual_amplification at this station.
function _amplification_residual(states, p)
    U2 = eltype(states).(p.U2)
    U2[3] = states[1]
    
    return [residual_amplification(p.xs, p.U1, U2, p.thermo)]
end

# 1-D Newton for amplification only.
function local_newton_amp(F, x0, p; tol=1e-10, max_iter=50)
    x = copy(x0)
    for i in 1:max_iter
        Fx = F(x, p)
        if norm(Fx) < tol
            return x, true
        end
        J = ForwardDiff.jacobian(s -> F(s, p), x)
        Δx = pinv(J) * (-Fx)
        # Damp step
        dm = 0.5 * (1.01 - i / max_iter)
        omega = abs(Δx[1]) > dm ? dm / abs(Δx[1]) : 1.0
        x .+= Δx .* omega
        if !p.turbulent && x[1] < 0
            x[1] = 0
        end
    end
    Fx = F(x, p)
    
    return x, norm(Fx) < tol
end

# ── Global Newton ───────────────────────────────────────────────────────────

"""
    global_newton(gr, U0; tol=1e-10, max_iter=50, verbose=true)

Coupled-viscous global Newton iteration. `gr` is built by
`build_global_residuals` and contains mutable arrays/refs that are
updated each iteration (stagnation index, transition index, turbulent
nodes, D matrix, ueinv signs, transition x's).

Returns `(U_final, converged)` where `U_final` is the flat state vector.
"""
function global_newton(gr, U0; tol=1e-10, max_iter=50, verbose=true)
    F_wrapper(x) = global_residual(x, gr)
    U = copy(vec(U0))

    num_nodes_total = length(gr.d)
    num_node_states = num_nodes_total * 4
    idx_thetas = 1:4:num_node_states
    idx_deltas = 2:4:num_node_states
    idx_state3s = 3:4:num_node_states
    idx_ues = 4:4:num_node_states

    for i in 1:max_iter
        R = F_wrapper(U)
        R_norm = norm(R[1:(3 * length(idx_ues))])

        if verbose
            println("Newton iteration $i: Rnorm = $R_norm")
        end

        if R_norm < tol
            if verbose
                println("Converged at step $i (R_norm = $R_norm)")
            end
            return U, true
        end

        deltas = U[idx_deltas]
        thetas = U[idx_thetas]
        state3s = U[idx_state3s]
        ues = U[idx_ues]
        ctmax = maximum(state3s[gr.turbulent_nodes])

        cj = ForwardDiff.JacobianConfig(F_wrapper, U, ForwardDiff.Chunk{4}())
        J = ForwardDiff.jacobian(F_wrapper, U, cj)
        ΔU = pinv(J) * (-R)

        # Zero out tiny state3 increments at nodes where state3 is exactly zero
        for j in idx_state3s
            if U[j] == 0 && abs(ΔU[j]) < 1e-8
                ΔU[j] = 0
            end
        end

        Δdeltas = ΔU[idx_deltas]
        Δthetas = ΔU[idx_thetas]
        Δstate3s = ΔU[idx_state3s]
        Δues = ΔU[idx_ues]
        omega = 1.0

        # Under-relaxation: cap fractional changes in each state group
        fmin = minimum(Δthetas ./ thetas)
        if fmin < -0.5
            omega = min(omega, abs(0.5 / fmin))
        end
        fmin = minimum(Δdeltas ./ deltas)
        if fmin < -0.5
            omega = min(omega, abs(0.5 / fmin))
        end

        # Prevent state3 from going negative on transitioned nodes
        for j in 1:length(state3s)
            !gr.turbulent_nodes[j] && state3s[j] < 0.2 && continue
            gr.turbulent_nodes[j] && state3s[j] < 0.1 * ctmax && continue
            (state3s[j] == 0 || Δstate3s[j] == 0) && continue
            if state3s[j] + Δstate3s[j] < 0
                omega = min(omega, 0.8 * abs(state3s[j] / Δstate3s[j]))
            end
        end

        # Cap changes in n / cτ / ue
        dnmax = maximum(abs.(Δstate3s[.!gr.turbulent_nodes]))
        if dnmax > 0
            omega = min(omega, abs(2 / dnmax))
        end
        dctmax = maximum(abs.(Δstate3s[gr.turbulent_nodes]))
        if dctmax > 0
            omega = min(omega, abs(0.5 / dctmax))
        end
        fmax = maximum(abs.(Δues) / gr.Vinf)
        if fmax > 0
            omega = min(omega, abs(0.2 / fmax))
        end

        U .+= ΔU .* omega

        # Repair bad Hk after update (push δ* up if Hk too low)
        for j in 1:length(thetas)
            Hkmin = j < gr.idx_w[1] ? 1.02 : 1.00005
            Hk = get_Hk(
                U[idx_deltas[j]],
                U[idx_thetas[j]],
                U[idx_ues[j]],
                gr.thermo;
                limit_Hk_max=false, limit_Hk_min=false,
            )
            if Hk < Hkmin
                U[idx_deltas[j]] += 2 * (Hkmin - Hk) * U[idx_thetas[j]]
            end
        end

        # Replace negative cτ on turbulent nodes with a small positive value
        for j in 1:length(state3s)
            if gr.turbulent_nodes[j] && U[idx_state3s[j]] < 0
                U[idx_state3s[j]] = 0.1 * ctmax
            end
        end

        # Stagnation move
        _, _, _, _, ues_new = stagnation_point_move!(gr, U[idx_ues])
        U[idx_ues] .= ues_new

        # Transition relocation
        Umat = reshape(U, 4, :)' |> Matrix
        resolve_transition_location!(gr, Umat; upper=false)
        resolve_transition_location!(gr, Umat; upper=true)
        U .= vec(Umat')
    end

    R = F_wrapper(U)
    R_norm = norm(R[1:(3 * length(idx_ues))])
    if verbose
        println(R_norm < tol ? "Converged at last step (R_norm = $R_norm)" :
                                "Failed to converge (R_norm = $R_norm)")
    end
    
    return U, R_norm < tol
end
