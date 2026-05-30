#=
Viscous-coupling orchestrator. Wraps the inviscid pipeline already
executed by `analyze` and adds wake construction, boundary-layer
initialization, the coupled global Newton, and final post-processing.

Single-body only; loops over the flow_angles vector externally.
=#

"""
    analyze_viscous(method::Xfoil, panel_geometry, system_geometry, system_matrices, strengths, flow_angles)

Entry point for the coupled viscous + inviscid Xfoil-style solve. Called
by `analyze` when `method.viscous == true`. Loops over each flow angle and
runs an independent viscous solve per α (single-body only).

# Returns
NamedTuple with both inviscid intermediates and per-α viscous results:
- Inviscid: `gamma_ref`, `A`.
- Per-α: `cl`, `cd`, `cm`, `cp`, `vs`, `converged`, `transition_upper`,
  `transition_lower`, plus per-surface BL state (`deltas_u/l/w`,
  `thetas_u/l/w`, `state3s_u/l/w`, `Veu/Vel/Vwake`).
"""
function analyze_viscous(
    method::Xfoil,
    panel_geometry,
    system_geometry,
    system_matrices,
    strengths,
    flow_angles,
)
    if system_geometry.nbodies != 1
        error("Viscous analysis currently supports a single body only; got $(system_geometry.nbodies).")
    end

    naoa = length(flow_angles)
    results = [
        _solve_one_alpha(method, panel_geometry, system_geometry, system_matrices,
                         strengths, flow_angles[a], a)
        for a in 1:naoa
    ]

    # Stack scalar results across alphas
    cl_arr = [r.cl for r in results]
    cd_arr = [r.cd for r in results]
    cm_arr = [r.cm for r in results]
    converged_arr = [r.converged for r in results]
    transition_upper_arr = [r.transition_upper for r in results]
    transition_lower_arr = [r.transition_lower for r in results]

    return (;
        # Inviscid intermediates
        gamma_ref=strengths[1:(panel_geometry.npanels + 1), :],
        A=system_matrices.A,
        # Per-α scalars
        cl=cl_arr,
        cd=cd_arr,
        cm=cm_arr,
        converged=converged_arr,
        transition_upper=transition_upper_arr,
        transition_lower=transition_lower_arr,
        # Per-α detailed results
        alphas=collect(flow_angles),
        per_alpha=results,
    )
end

# Run the full coupled viscous solve for a single angle of attack `alpha_deg`.
function _solve_one_alpha(
    method, panel_geometry, system_geometry, system_matrices, strengths, alpha_deg, alpha_idx
)
    alpha = alpha_deg * π / 180   # FLOWFoil API uses degrees externally
    cosalpha = cos(alpha)
    sinalpha = sin(alpha)
    chord = method.chord
    Vinf = method.Vinf
    nnodes = panel_geometry.npanels + 1
    gamma_ref = strengths[1:nnodes, :]

    # ── Build viscous-side geometry ─────────────────────────────────────────
    airfoil_panels = build_viscous_airfoil_panels(panel_geometry)
    te_geometry = build_te_geometry(airfoil_panels)

    # Tangential velocity for find_stagnation_split (signed) and ue_inv (absolute on airfoil)
    Vt_signed = gamma_ref * [cosalpha, sinalpha]
    ps = find_stagnation_split(
        panel_geometry.nodes[:, 1],
        panel_geometry.nodes[:, 2],
        Vt_signed,
        airfoil_panels,
        method.xft_l,
        method.xft_u,
    )

    numpanels_wake =
        method.numpanels_wake === nothing ? Int(ceil(nnodes / 10 + 10) - 1) :
        method.numpanels_wake

    wake_panels, ue_wake_ref, wake_gap = build_wake(
        airfoil_panels,
        te_geometry,
        ps.s,
        gamma_ref,
        alpha,
        Vinf,
        chord,
        method.wake_length,
        numpanels_wake,
    )

    # Combined inviscid edge velocity at every node (airfoil: |Vt|, wake: signed)
    Vt_abs = abs.(Vt_signed)
    Vwake_inv = ue_wake_ref * [cosalpha, sinalpha]
    Vwake_inv[1] = Vt_abs[end]   # match airfoil last node
    ueinv = vcat(Vt_abs, Vwake_inv)

    # Arc-length vectors
    _, sw, swake0 = arclength_vectors(airfoil_panels, wake_panels, ps.s_stag)
    swake0 = swake0 .- swake0[1] .+ ps.s[end]
    idx_w = collect((ps.idx_u[end] + 1):(ps.idx_u[end] + length(sw)))

    # ── Influence matrices ──────────────────────────────────────────────────
    D, d, Bprime, Dwake = create_D_influence_matrix(
        airfoil_panels, te_geometry, wake_panels, system_matrices.A, Vt_signed
    )

    # ── Thermo + state refs ─────────────────────────────────────────────────
    thermo = init_viscous_thermo(
        method.Re, method.Mach, method.rho, chord, Vinf; ncrit=method.ncrit
    )
    state_refs = (;
        sft_l=Ref(ps.sft_l),
        sft_u=Ref(ps.sft_u),
        transition_lower=Ref(chord),
        transition_upper=Ref(chord),
        wake_gap=wake_gap,
    )

    # ── BL initialization on each surface ──────────────────────────────────
    if method.verbose
        println("\nInitializing boundary layer (α = $(alpha_deg) deg)...")
    end

    Vel = ueinv[ps.idx_l[1]:-1:ps.idx_l[end]]
    Veu = ueinv[ps.idx_u[1]:ps.idx_u[end]]
    Vew = ueinv[idx_w]

    bl_l = initialize_boundary_layer(
        ps.sl, thermo, state_refs, copy(Vel); upper=false, airfoil=true,
        max_iter=method.maxiters, tol=method.etol,
    )
    bl_u = initialize_boundary_layer(
        ps.su, thermo, state_refs, copy(Veu); upper=true, airfoil=true,
        max_iter=method.maxiters, tol=method.etol,
    )

    is_turb_l = bl_l.idx_transition > length(bl_l.deltas) ?
                falses(length(ps.sl)) :
                [i >= bl_l.idx_transition for i in 1:length(ps.sl)]
    is_turb_u = bl_u.idx_transition > length(bl_u.deltas) ?
                falses(length(ps.su)) :
                [i >= bl_u.idx_transition for i in 1:length(ps.su)]

    # Assemble initial global state matrix
    num_states_total = nnodes + length(sw)
    U_global_mat = zeros(num_states_total, 4)
    U_global_mat[ps.idx_l[1]:-1:ps.idx_l[end], 1] .= bl_l.thetas
    U_global_mat[ps.idx_l[1]:-1:ps.idx_l[end], 2] .= bl_l.deltas
    U_global_mat[ps.idx_l[1]:-1:ps.idx_l[end], 3] .= bl_l.state3s
    U_global_mat[ps.idx_l[1]:-1:ps.idx_l[end], 4] .= bl_l.ues
    U_global_mat[ps.idx_u[1]:ps.idx_u[end], 1] .= bl_u.thetas
    U_global_mat[ps.idx_u[1]:ps.idx_u[end], 2] .= bl_u.deltas
    U_global_mat[ps.idx_u[1]:ps.idx_u[end], 3] .= bl_u.state3s
    U_global_mat[ps.idx_u[1]:ps.idx_u[end], 4] .= bl_u.ues

    # Wake initialization (first node bridges from TE)
    U_stag_wake = initialize_wake_first_node(
        U_global_mat[ps.idx_l[end], :],
        U_global_mat[ps.idx_u[end], :],
        zeros(3),
        te_geometry,
        is_turb_l[end],
        is_turb_u[end],
        thermo,
        state_refs;
        init_flag=true,
    )
    bl_w = initialize_boundary_layer(
        sw, thermo, state_refs, copy(Vew), U_stag_wake;
        airfoil=false, max_iter=method.maxiters, tol=method.etol,
    )
    U_global_mat[idx_w, 1] .= bl_w.thetas
    U_global_mat[idx_w, 2] .= bl_w.deltas
    U_global_mat[idx_w, 3] .= bl_w.state3s
    U_global_mat[idx_w, 4] .= bl_w.ues

    # Per-node turbulent flag
    turbulent_nodes = fill(false, num_states_total)
    turbulent_nodes[ps.idx_l[1]:-1:ps.idx_l[end]] .= is_turb_l
    turbulent_nodes[ps.idx_u[1]:ps.idx_u[end]] .= is_turb_u
    turbulent_nodes[idx_w] .= true

    # ── Build GlobalResiduals and solve ────────────────────────────────────
    gr = build_global_residuals(;
        s=ps.s, su=ps.su, sl=ps.sl, sw=sw,
        idx_u=[ps.idx_u[1], ps.idx_u[end]],
        idx_l=[ps.idx_l[1], ps.idx_l[end]],
        idx_w=idx_w,
        idx_transition_u=[bl_u.idx_transition],
        idx_transition_l=[bl_l.idx_transition],
        turbulent_nodes=turbulent_nodes,
        thermo=thermo, state_refs=state_refs,
        Bprime=Bprime, Dwake=Dwake, D=D, d=d, ueinv=ueinv,
        airfoil_panels=airfoil_panels, te_geometry=te_geometry,
        wake_panels=wake_panels, swake0=swake0,
        Vinf=Vinf, Mach=method.Mach,
        verbose=method.verbose,
    )

    U0 = vec(U_global_mat')
    solution, converged = global_newton(
        gr, U0; tol=method.etol, max_iter=method.maxiters, verbose=method.verbose
    )

    if method.verbose && !converged
        @warn "Global viscous solve did not converge (α = $(alpha_deg) deg)."
    end

    U_final_mat = reshape(solution, 4, :)' |> Matrix

    # ── Post-process ─────────────────────────────────────────────────────────
    thetas_final = U_final_mat[:, 1]
    deltas_final = U_final_mat[:, 2]
    state3s_final = U_final_mat[:, 3]
    ues_final = get_uk(U_final_mat[:, 4], thermo)

    # Surface slices
    deltas_u = deltas_final[ps.idx_u[1]:ps.idx_u[end]]
    deltas_l = deltas_final[ps.idx_l[1]:-1:ps.idx_l[end]]
    deltas_w = deltas_final[idx_w[1]:idx_w[end]]
    thetas_u = thetas_final[ps.idx_u[1]:ps.idx_u[end]]
    thetas_l = thetas_final[ps.idx_l[1]:-1:ps.idx_l[end]]
    thetas_w = thetas_final[idx_w[1]:idx_w[end]]
    state3s_u = state3s_final[ps.idx_u[1]:ps.idx_u[end]]
    state3s_l = state3s_final[ps.idx_l[1]:-1:ps.idx_l[end]]
    state3s_w = state3s_final[idx_w[1]:idx_w[end]]
    Veu = ues_final[ps.idx_u[1]:ps.idx_u[end]]
    Vel = ues_final[ps.idx_l[1]:-1:ps.idx_l[end]]
    Vwake = ues_final[idx_w[1]:idx_w[end]]

    # Cp (compressible Karman-Tsien correction if needed)
    cp = 1 .- (U_final_mat[:, 4] ./ Vinf) .^ 2
    if method.Mach > 0
        den = thermo.KTb .+ 0.5 .* thermo.KTl .* (1 + thermo.KTb) .* cp
        cp ./= den
    end

    # Integrated forces
    cl = _viscous_cl(panel_geometry, ues_final, ps.idx_u, Vinf, sinalpha, cosalpha)
    cd = _viscous_cd(thetas_w, Vwake, deltas_w, Vinf)
    cm = _viscous_cm(panel_geometry, system_geometry, cp, cosalpha, sinalpha)

    return (;
        cl=cl,
        cd=cd,
        cm=cm,
        cp=cp[1:nnodes],
        vs=U_final_mat[1:nnodes, 4],
        converged=converged,
        transition_upper=state_refs.transition_upper[],
        transition_lower=state_refs.transition_lower[],
        deltas_u=deltas_u, deltas_l=deltas_l, deltas_w=deltas_w,
        thetas_u=thetas_u, thetas_l=thetas_l, thetas_w=thetas_w,
        state3s_u=state3s_u, state3s_l=state3s_l, state3s_w=state3s_w,
        Veu=Veu, Vel=Vel, Vwake=Vwake,
        s=ps.s, su=ps.su, sl=ps.sl, sw=sw,
        idx_u=ps.idx_u, idx_l=ps.idx_l, idx_w=idx_w,
    )
end

# Squire-Young-style cl from the airfoil edge velocities
function _viscous_cl(panel_geometry, ues_final, idx_u, Vinf, sinalpha, cosalpha)
    nnodes = size(panel_geometry.nodes, 1)
    xs = panel_geometry.nodes[:, 1]
    ys = panel_geometry.nodes[:, 2]
    ues = ues_final[1:idx_u[2]]
    cps = 1 .- (ues ./ Vinf) .^ 2

    cl = zero(eltype(cps))
    cn = zero(eltype(cps))
    ct = zero(eltype(cps))
    for i in 1:(length(ues) - 1)
        Δx = xs[i + 1] - xs[i]
        Δy = ys[i + 1] - ys[i]
        cpbar = (cps[i] + cps[i + 1]) / 2
        cn -= cpbar * Δx     # normal force (per chord)
        ct += cpbar * Δy     # tangential force
    end

    return cn * cosalpha - ct * sinalpha
end

# Squire-Young-style cd from the wake far-field state.
function _viscous_cd(thetas_w, Vwake, deltas_w, Vinf)
    H = deltas_w[end] / thetas_w[end]
    ue_te = Vwake[end]

    # Squire-Young: cd = 2 · θ_∞ · (ue_∞ / V_∞)^((5+H)/2)
    return 2 * thetas_w[end] * (ue_te / Vinf)^((5 + H) / 2)
end

# Pitching moment coefficient at quarter chord from the Cp distribution.
function _viscous_cm(panel_geometry, system_geometry, cp, cosalpha, sinalpha)
    nnodes = panel_geometry.npanels + 1
    cp_af = cp[1:nnodes]
    chord = system_geometry.chord_length
    x0 = chord / 4
    z0 = 0.0
    panel_vector = panel_geometry.panel_vectors
    panel_edges = panel_geometry.panel_edges
    cmmat = [2 1; 1 2] ./ 6

    dxddmi = [
        panel_vector[i, 1] * (panel_edges[i, 1, 1] - x0) +
        panel_vector[i, 2] * (panel_edges[i, 1, 2] - z0)
        for i in 1:panel_geometry.npanels
    ]
    dxddmip1 = [
        panel_vector[i, 1] * (panel_edges[i, 2, 1] - x0) +
        panel_vector[i, 2] * (panel_edges[i, 2, 2] - z0)
        for i in 1:panel_geometry.npanels
    ]
    
    return sum(
        ([cp_af[i] cp_af[i + 1]] * cmmat * [dxddmi[i]; dxddmip1[i]])[1]
        for i in 1:panel_geometry.npanels
    ) / chord^2
end
