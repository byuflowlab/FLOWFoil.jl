function run_coupled_system_global(xaf::AbstractVector, yaf::AbstractVector; 
                                    viscous::Bool=true, alpha=0.0, Re=1e6, Mach=0.0, 
                                    rho=1.0, chord=1.0, Vinf=1.0, inviscidproblem=XFOIL()::InviscidProblem, 
                                    ncrit=9.0,
                                    etol=1e-10,
                                    maxiters=50,
                                    name="", 
                                    verbose=true, plot_bool=true, 
                                    cosine_spacing_flag=false, 
                                    data_file=nothing, # for plotting comparison
                                    use_native_solver_local=true,
                                    use_native_solver_global=true,
                                    hard_code_geometry=false,
                                    hard_code_initial_states=false,
                                    transition_method=1,
                                    method=Xfoil(),
                                    kwargs...)

    coordinates = [xaf yaf]
    flow_angles = alpha
    coordinates, nbodies, flow_angles = reformat_inputs(coordinates, flow_angles)

    if nbodies > 1
        println("Viscous can only work with one body for now")
    end

    panel_geometry = generate_panel_geometry(method, coordinates)
    """
    - Named tuple `panel_geometry` containing:
    - `npanels::Int`: Number of panels (nodes - 1).
    - `panel_edges::Array{TF, 3}`: An (npanels, 2, 2) array containing the start and end points of each panel.
    - `panel_vectors::Array{TF, 2}`: An (npanels, 2) array of vectors describing panel directions.
    - `panel_lengths::Array{TF}`: Vector of panel lengths.
    - `nodes::Array{TF, 2}`: Array of node coordinates (npanels + 1, 2).
    """
    system_geometry = generate_system_geometry(method, panel_geometry)

    """
    - A named tuple `system_geometry` containing:
    - `nbodies`: Number of bodies/airfoils.
    - `panel_indices`: Vector of panel index ranges for each body.
    - `node_indices`: Vector of node index ranges for each body.
    - `nodes`: Combined matrix of all nodes for all bodies.
    - `mesh2panel`: Mapping from global mesh index to panel index.
    - `chord_length`: Overall chord length of the system (max trailing edge - min leading edge).
    - `panel_length`: Vector of panel lengths for all panels combined.
    - `r1`, `lnr1`, `r1normal`, `r1tangent`, `theta1`, `r2`, `lnr2`, `theta2`: Matrices of precomputed influence geometry quantities for all panels on all field points.
    - `TE_geometry`: Named tuple with trailing edge gap panel geometry and related influence parameters:
        - `blunt_te`: Boolean vector indicating blunt trailing edges.
        - `panel_length`: Trailing edge gap lengths.
        - `tdp`, `txp`: Trailing edge panel parameters.
        - `r1`, `lnr1`, `r1normal`, `r1tangent`, `theta1`, `r2`, `lnr2`, `theta2`: Influence geometry matrices for trailing edge panels.
    """

    airfoil = Airfoil("", xaf, yaf, system_geometry.TE_geometry.blunt_te[1])

    mu, KTb, KTl, H0, mu0, rho0 = initialize_thermo(Re, Mach, rho, chord, Vinf)

    operatingparameters = OperatingParameters(; numpanels_airfoil=panel_geometry.npanels, blunt_TE=system_geometry.TE_geometry.blunt_te[1], 
                                                Re, Mach, rho, mu, Vinf, chord, alpha=flow_angles, 
                                                viscous, ncrit, etol, maxiters,
                                                H0, mu0, rho0, KTb, KTl,
                                                use_native_solver_local, verbose)

    panels_airfoil = get_panels(airfoil, panel_geometry.npanels; cosine_spacing_flag)

    if verbose; println("Running Inviscid solution . . ."); end
    system_matrices = generate_system_matrices(method, panel_geometry, system_geometry)
    system_matrices = (; A=system_matrices.A, b=system_matrices.b * Vinf)
    """
    - Named tuple with:
    - `A`: The coefficient matrix (influence matrix) representing panel interactions.
    - `b`: The boundary conditions vector representing the RHS of the linear system.
    - `node_indices`: The node indexing ranges for each body in the system.
    """

    strengths = solve(method, system_matrices)
    """
     - `x::inviscid solution, consisting of node vortex strengths`
    """

    post_processed = post_process(method, panel_geometry, system_geometry, strengths, flow_angles)

    β = zeros(panel_geometry.npanels+1,panel_geometry.npanels)
    xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[end].x2)
    ys = vcat([p.y1 for p in panels_airfoil], panels_airfoil[end].y2)
    for i in 1:panel_geometry.npanels+1
        for j in 1:panel_geometry.npanels
            β[i,j] = get_beta_mfoil(i, j, xs[i], ys[i], panels_airfoil[j])
        end
    end

    inviscidoutputs = (; cl=vec(post_processed.cl), cl_KJ=0, cd=vec(post_processed.cd), cmle=0,
                        cmqc=vec(post_processed.cm), cps=post_processed.cp, A=vcat(-1 .* system_matrices.A[1:end-2, :], system_matrices.A[end-1:end, :]),
                        b=-1 .* system_matrices.b, γref=strengths[1:end-1, :], r=hcat(system_geometry.r1, system_geometry.r2[:, end]), 
                        β= β, xstar=system_geometry.r1tangent,
                        ystar=system_geometry.r1normal, logij=hcat(system_geometry.lnr1, system_geometry.lnr2[:, end]),
                        solution=post_processed.vs, velocities=post_processed.vs)
    #cd is different. they don't have clKJ, cmle. beta should be probably caclulated as theta2-theta1. The diagonal and corner entries were not working with that
    #(; cl, cl_KJ, cd, cmle, cmqc, cps, A, b, γref, r, β, xstar, ystar, logij, solution, velocities)

    # return inviscidoutputs
    if !viscous
        viscousoutputs = nothing
        return (; airfoil, panels_airfoil, operatingparameters, inviscidoutputs, viscousoutputs) 
    end

    viscousoutputs = run_viscous(airfoil, panels_airfoil, operatingparameters, inviscidoutputs, sin(alpha), cos(alpha); verbose)

    if plot_bool
        which_plots = [:geometry, :cp, :cp_comparison, :blayers, :blayers_paper, :cp_paper]
        plotting_functions(panels_airfoil, inviscidoutputs, viscousoutputs, which_plots; verbose, alphas=operatingparameters.alpha, Vinfs=operatingparameters.Vinf, data_file)
    end

    return (; airfoil, panels_airfoil, operatingparameters, inviscidoutputs, viscousoutputs)  
end

function run_viscous(airfoil::Airfoil,
                    panels_airfoil,
                    operatingparameters::OperatingParameters,
                    inviscidoutputs, 
                    sinalpha, cosalpha; 
                    name="", 
                    verbose=true, 
                    hard_code_geometry=false,
                    hard_code_initial_states=false,
                    transition_method=1,
                    reinit=true,
                    initialized_states=nothing,
                    detailed_output_flag=true,
                    use_native_solver_global=true,
                    alpha_idx=1,
                    kwargs...)

    panel_TE = get_TE_panel(panels_airfoil)

    Vairfoil = inviscid_velocities_airfoil(sinalpha, cosalpha, inviscidoutputs.γref)

    # get s, ue for each surface, stagnation point
    #todo: check to make sure airfoil.x is the same thing as xaf - do I change one based on cosine spacing? If I do that, do I update airfoil.x and airfoil.y?
    s, su, idx_u, sl, idx_l, s_stag = parseairfoil(airfoil.x, airfoil.y, Vairfoil, panels_airfoil, operatingparameters; verbose=false) #assume single alpha value for now
    
    # build wake
    panels_wake, ue_wake_ref = build_wake(operatingparameters, panels_airfoil, panel_TE, s, inviscidoutputs.γref, sinalpha, cosalpha) #update to make the wake follow the streamline

    ue_inv_solution = inviscid_velocities_airfoil_and_wake(sinalpha, cosalpha, inviscidoutputs.γref, ue_wake_ref)

    # wake s, idx vectors
    _, sw, swake0 = get_s_length_vectors(panels_airfoil, panels_wake, s_stag)
    swake0 = swake0 .- swake0[1] .+ s[end]
    idx_w = range(idx_u[end]+1, step=1, length=length(sw))
    
    # edge velocities
    Vel = ue_inv_solution[idx_l[1]:-1:idx_l[end]]
    Veu = ue_inv_solution[idx_u[1]:idx_u[end]]
    Vew = ue_inv_solution[idx_w]

    D, d, Bprime, Dwake = create_D_influence_matrix(panels_airfoil, panel_TE, panels_wake, inviscidoutputs; alpha_idx)

    if reinit
        # initialize boundary layer on upper and lower surfaces
        if verbose; println("\nInitializing boundary layer surfaces . . ."); end
        deltas_l, thetas_l, state3s_l, Vel, 
            idx_transition_l, xt_l, theta_t_l, delta_t_l = initialize_boundary_layer(sl, operatingparameters, Vel; upper=false, airfoil=true, verbose)
        deltas_u, thetas_u, state3s_u, Veu, 
            idx_transition_u, xt_u, theta_t_u, delta_t_u = initialize_boundary_layer(su, operatingparameters, Veu; upper=true, airfoil=true, verbose)

        is_turbulent_l = idx_transition_l == length(deltas_l)+1 ? zeros(Bool, length(sl)) : [i >= idx_transition_l for i in 1:length(sl)]
        is_turbulent_u = idx_transition_u == length(deltas_u)+1 ? zeros(Bool, length(su)) : [i >= idx_transition_u for i in 1:length(su)]

        # initialize global residuals
        num_states = operatingparameters.numpanels_airfoil+1+operatingparameters.numpanels_wake+1 #total # of nodes
            
        U_global_0_temp = zeros(num_states, 4) 

        # assign to global state array
        U_global_0_temp[idx_l[1]:-1:idx_l[end],1] = thetas_l
        U_global_0_temp[idx_l[1]:-1:idx_l[end],2] = deltas_l
        U_global_0_temp[idx_l[1]:-1:idx_l[end],3] = state3s_l
        U_global_0_temp[idx_l[1]:-1:idx_l[end],4] = Vel

        U_global_0_temp[idx_u[1]:idx_u[end],1] = thetas_u
        U_global_0_temp[idx_u[1]:idx_u[end],2] = deltas_u
        U_global_0_temp[idx_u[1]:idx_u[end],3] = state3s_u
        U_global_0_temp[idx_u[1]:idx_u[end],4] = Veu

        # # initialize boundary layer on wake
        U_stag_wake = initialize_wake_first_node(U_global_0_temp[idx_l[end],:], U_global_0_temp[idx_u[end],:], zeros(3),
                                                                panel_TE, is_turbulent_l[end], is_turbulent_u[end],
                                                                operatingparameters; init_flag=true)

        deltas_w, thetas_w, state3s_w, Vew, _, _, _, _ = initialize_boundary_layer(sw, operatingparameters, Vew, U_stag_wake; airfoil=false, verbose)

        is_turbulent_l = idx_transition_l == length(deltas_l)+1 ? zeros(Bool, length(sl)) : [i >= idx_transition_l for i in 1:length(sl)]
        is_turbulent_u = idx_transition_u == length(deltas_u)+1 ? zeros(Bool, length(su)) : [i >= idx_transition_u for i in 1:length(su)]

        # assign wake states to global array
        U_global_0_temp[idx_w,1] = thetas_w
        U_global_0_temp[idx_w,2] = deltas_w
        U_global_0_temp[idx_w,3] = state3s_w
        U_global_0_temp[idx_w,4] = Vew

        turbulent_nodes = fill(false, idx_w[end])
        turbulent_nodes[idx_l[1]:-1:idx_l[end]] .= is_turbulent_l
        turbulent_nodes[idx_u[1]:idx_u[end]] .= is_turbulent_u
        turbulent_nodes[idx_w] .= true

        #state array to state vector
        U_global_0_temp = reshape(U_global_0_temp', :, 1)

        if operatingparameters.transition_method == 3
            # assign xt, theta_t, and delta_t states
            U_global_0_transition = vcat(xt_l, theta_t_l, delta_t_l,
                                        xt_u, theta_t_u, delta_t_u)
            U_global_0 = vcat(U_global_0_temp, U_global_0_transition)
        else
            U_global_0 = U_global_0_temp
        end
    else
        @assert !isnothing(initialized_states) "if reinit is false, initialized states must be provided."
        U_global_0 = initialized_states.U
        idx_transition_l = initialized_states.idx_transition_l
        idx_transition_u = initialized_states.idx_transition_u
        turbulent_nodes = initialized_states.turbulent_nodes
        idx_l = initialized_states.idx_l
        idx_u = initialized_states.idx_u
        idx_w = initialized_states.idx_w
    end

    global_residuals = GlobalResiduals(s, sl, su, sw, idx_u, idx_l, idx_w, 
                                    [idx_transition_u], [idx_transition_l], 
                                    turbulent_nodes, operatingparameters, Bprime, Dwake, D, d, 
                                    ue_inv_solution,
                                    panels_airfoil, panel_TE, panels_wake,
                                    inviscidoutputs, swake0)

    if !use_native_solver_global
        prob = NonlinearSolve.NonlinearProblem(global_residuals, U_global_0)

        isoutofdomain = IsOutOfDomain(fill(1e-6, 3)) # prevent Newton solve from going to negative values #todo: generalize better or make a user input?
        sol = NonlinearSolve.solve(prob, NonlinearSolve.LevenbergMarquardt(), 
                                        abstol = operatingparameters.etol, 
                                        maxiters = operatingparameters.maxiters,
                                        isoutofdomain = isoutofdomain)
                                        
        if sol.retcode != NonlinearSolve.ReturnCode.Success
            prob2 = NonlinearSolve.NonlinearProblem(global_residuals, U_global_0)
            
            sol2 = NonlinearSolve.solve(prob2, NonlinearSolve.NewtonRaphson(), 
                                        abstol = operatingparameters.etol, 
                                        maxiters = operatingparameters.maxiters,
                                        isoutofdomain = isoutofdomain)

            if sol2.retcode == NonlinearSolve.ReturnCode.Success
                sol = sol2
            end
        end
        if verbose; println("Global Newton convergence: ", sol.retcode); end
        U_global_final = reshape(sol.u, 4, :)'
    else
        # solution, converged = global_levenberg_marquardt(global_residuals, U_global_0; tol=etol, max_iter=maxiters)
        solution, converged = global_newton(global_residuals, U_global_0; 
                                            tol=operatingparameters.etol, 
                                            max_iter=operatingparameters.maxiters, 
                                            verbose)
        if operatingparameters.transition_method == 3
            U_global_final = reshape(solution[1:end-6], 4, :)'
        else
            U_global_final = reshape(solution, 4, :)'
        end
    end

    if !converged && verbose
        @warn "Global solve did not converge." 
    end

    thetas_final = U_global_final[:,1]
    deltas_final = U_global_final[:,2]
    state3s_final = U_global_final[:,3]
    ues_final = get_uk(U_global_final[:,4], operatingparameters)

    deltas_u = deltas_final[idx_u[1]:idx_u[end]]
    deltas_l = deltas_final[idx_l[1]:-1:idx_l[end]]
    deltas_w = deltas_final[idx_w[1]:idx_w[end]]
    thetas_u = thetas_final[idx_u[1]:idx_u[end]]
    thetas_l = thetas_final[idx_l[1]:-1:idx_l[end]]
    thetas_w = thetas_final[idx_w[1]:idx_w[end]]
    state3s_u = state3s_final[idx_u[1]:idx_u[end]]
    state3s_l = state3s_final[idx_l[1]:-1:idx_l[end]]
    state3s_w = state3s_final[idx_w[1]:idx_w[end]]
    Veu = ues_final[idx_u[1]:idx_u[end]]
    Vel = ues_final[idx_l[1]:-1:idx_l[end]]
    Vwake = ues_final[idx_w[1]:idx_w[end]]

    su, sl, sw, s_stag, ues_final = stagnation_point_move!(global_residuals, ues_final)

    cps = 1 .- (U_global_final[:,4] ./ operatingparameters.Vinf) .^ 2
    if operatingparameters.Mach > 0
        den = operatingparameters.KTb .+ 0.5 .* operatingparameters.KTl .* (1+operatingparameters.KTb) .* cps
        cps = cps./den
    end

    stu = operatingparameters.transition_upper[]
    stl = operatingparameters.transition_lower[]

    idx_transition_l = global_residuals.idx_transition_l[1]
    idx_transition_u = global_residuals.idx_transition_u[1]

    transition_upper = operatingparameters.transition_upper
    transition_lower = operatingparameters.transition_lower

    viscousoutputs = (; deltas_u, deltas_l, deltas_w,
                thetas_u, thetas_l, thetas_w,
                state3s_u, state3s_l, state3s_w,
                Veu, Vel, Vwake, cps, s, su, sl, sw,
                idx_u, idx_l, idx_w, stu, stl, s_stag, 
                idx_transition_l, idx_transition_u)

    # Integrated outputs from Anna
    cl = get_cl(airfoil.x, airfoil.y, ues_final[1:idx_u[2]], operatingparameters.Vinf, sinalpha, cosalpha)
    cd = get_cd(viscousoutputs, operatingparameters)
    # cdp, cd = get_cdp(panels_airfoil, viscousoutputs, operatingparameters)

    if detailed_output_flag
        viscousoutputs = (; deltas_u, deltas_l, deltas_w,
                thetas_u, thetas_l, thetas_w,
                state3s_u, state3s_l, state3s_w,
                Veu, Vel, Vwake, cps, s, su, sl, sw,
                idx_u, idx_l, idx_w, stu, stl, s_stag, 
                idx_transition_l, idx_transition_u,
                cl, cd,
                transition_upper, transition_lower,
                panels_wake, panel_TE, converged)
        return viscousoutputs
    else
        U_global_final = reshape(U_global_final', :, 1)[:]
        BLinitializedstates = BoundaryLayerInitializedStates(U_global_final, idx_transition_l, idx_transition_u, turbulent_nodes, idx_l, idx_u, idx_w)
        return (; cl, cd, converged, BLinitializedstates)
    end
end

struct GlobalResiduals{TF}
    s::AbstractVector{TF}
    sl::AbstractVector{TF}
    su::AbstractVector{TF}
    sw::AbstractVector{TF}
    idx_u::AbstractVector{Int64}
    idx_l::AbstractVector{Int64}
    idx_w::AbstractVector{Int64}
    idx_transition_u::AbstractVector{Int64}
    idx_transition_l::AbstractVector{Int64}
    turbulent_nodes::AbstractVector{Bool}
    operatingparameters::OperatingParameters
    Bprime::AbstractArray{TF}
    Dwake::AbstractArray{TF}
    D::AbstractArray{TF}
    d::AbstractVector{TF}
    ueinv::AbstractVector{TF}
    panels_airfoil::AbstractVector{Panel{TF}}
    panel_TE::TEgeometry{TF}
    panels_wake::AbstractVector{WakePanel{TF}}
    inviscidoutputs::NamedTuple #todo: check on this
    swake0::AbstractVector{TF}
end

function (p::GlobalResiduals)(U, params)

    (; s, idx_u, idx_l, idx_w, operatingparameters, D, ueinv, 
        idx_transition_u, idx_transition_l, d, 
        panels_wake, panel_TE, turbulent_nodes, swake0) = p

    if operatingparameters.transition_method == 3
        U_nodes = reshape(U[1:end-6], 4, :)'
        U_transition = U[end-5:end]
    else
        U_nodes = reshape(U, 4, :)'
    end
    Rs_global = zeros(eltype(U_nodes), length(ueinv), 3)
    
    if p.operatingparameters.transition_method == 3
        Rs_transition = zeros(eltype(U_nodes), 6)
    end

    # recalculate exact stagnation point
    su, sl, sw, _ = split_s(s, U_nodes[:,4], idx_u, idx_l, idx_w, swake0)

    ####################### Lower surface #######################
    upper = false
    turbulent = false
    transition = false
    airfoil = true
    ## check if stagnation is exactly on node
    if sl[1] < 1e-8*sl[end]
        i0 = 2 # hit the stagnation state exactly
        idx_stagnation = idx_l[1] - 1
    else
        i0 = 1
        idx_stagnation = idx_l[1]
    end

    # stagnation point treatment
        # extrapolate stagnation based on first two nodes in the boundary layer
    Ust, xst = interpolate_stagnation_state(U_nodes[idx_stagnation,:], U_nodes[idx_stagnation-1,:], sl[i0], sl[i0+1]) # method found at end of viscous file

    similar = true
    forced_transition_flag = false
    stagnation_residuals = ResidualParamsStation([xst, xst], Ust[2], Ust[1], Ust[3], [Ust[4], Ust[4]], 
                                            zeros(eltype(U),2), operatingparameters, upper,
                                            turbulent, transition, similar, airfoil, forced_transition_flag)

    Rs_global[idx_stagnation,:] .= stagnation_residuals(Ust[1:3], nothing)
    similar = false

    if i0 == 2
        Rs_global[idx_l[1],:] = U_nodes[idx_l[1],1:3] .- Ust[1:3]
    end

    idx_surface = idx_l[1]:-1:idx_l[end]

    for i in i0+1:length(idx_surface)
        idx = idx_surface[i]
           
        if !turbulent && (i == idx_transition_l[1] || sl[i] > operatingparameters.sft_l[1]) # for now, assume transition panel is the same as before
            transition = true
            if sl[i] > operatingparameters.sft_l[1] # forced transition
                forced_transition_flag = true
            end
        end

        station_residuals = ResidualParamsStation(eltype(U).(sl[i-1:i]), U_nodes[idx+1,2], U_nodes[idx+1,1], U_nodes[idx+1,3], U_nodes[(idx+1:-1:idx),4], 
                                                zeros(eltype(U),2), operatingparameters, upper, 
                                                turbulent, transition, similar, airfoil, forced_transition_flag)
        
        U_temp = (transition && operatingparameters.transition_method == 3) ? vcat(U_nodes[idx,1:3], U_transition[1:3]) : U_nodes[idx,1:3]
        Rs_local = station_residuals(U_temp, nothing)

        if !transition || p.operatingparameters.transition_method != 3
            Rs_global[idx,:] .= Rs_local
        # else
        #     Rs_global[idx,:] .= Rs_local[4:6]
        #     Rs_transition[1:3] .= Rs_local[1:3]
        end

        if transition # if we just did a transition station, turn it off and turn on turbulent
            transition = false
            forced_transition_flag = false
            turbulent = true
        end

        if turbulent
            turbulent_nodes[idx] = true
        end
    end

    TE_l_isturbulent = turbulent

    ####################### Upper surface #######################
    upper = true
    turbulent = false
    transition = false
    airfoil = true

    ## check if stagnation is exactly on node
    if su[1] < 1e-8*su[end]
        i0 = 2 # hit the stagnation state exactly
        idx_stagnation = idx_u[1] + 1
    else
        i0 = 1
        idx_stagnation = idx_u[1]
    end

    # stagnation point treatment - extrapolate stagnation based on first two nodes in the boundary layer
    Ust, xst = interpolate_stagnation_state(U_nodes[idx_stagnation,:], U_nodes[idx_stagnation+1,:], su[i0], su[i0+1]) # method found at end of viscous file

    similar = true
    forced_transition_flag = false
    stagnation_residuals = ResidualParamsStation([xst, xst], Ust[2], Ust[1], Ust[3], [Ust[4], Ust[4]], 
                                            zeros(eltype(U),2), operatingparameters, upper,
                                            turbulent, transition, similar, airfoil, forced_transition_flag)

    Rs_global[idx_stagnation,:] .= stagnation_residuals(Ust[1:3], nothing)
    similar = false

    if i0 == 2
        Rs_global[idx_u[1],:] = U_nodes[idx_u[1],1:3] .- Ust[1:3]
    end
    idx_surface = idx_u[1]:1:idx_u[end]

    for i in i0+1:length(idx_surface)
        idx = idx_surface[i]

        if !turbulent && (i == idx_transition_u[1] || su[i] > operatingparameters.sft_u[1]) # for now, assume transition panel is the same as before
            transition = true
            if su[i] > operatingparameters.sft_u[1]
                forced_transition_flag = true
            end
        end

        station_residuals = ResidualParamsStation(eltype(U).(su[i-1:i]), U_nodes[idx-1,2], U_nodes[idx-1,1], U_nodes[idx-1,3], U_nodes[(idx-1):idx,4], 
                                                    zeros(eltype(U),2), operatingparameters, upper, 
                                                    turbulent, transition, similar, airfoil, forced_transition_flag)

        U_temp = (transition && operatingparameters.transition_method == 3) ? vcat(U_nodes[idx,1:3], U_transition[4:6]) : U_nodes[idx,1:3]
        Rs_local = station_residuals(U_temp, nothing)
        
        if !transition || p.operatingparameters.transition_method != 3
            Rs_global[idx,:] .= Rs_local
        # else
        #     Rs_global[idx,:] .= Rs_local[4:6]
        #     Rs_transition[4:6] .= Rs_local[1:3]
        end

        if transition # if we just did a transition station, turn it off and turn on turbulent
            transition = false
            turbulent = true
            debug = false
            forced_transition_flag = false
        end

        if turbulent
            turbulent_nodes[idx] = true
        end
    end

    TE_u_isturbulent = turbulent

    turbulent = true
    airfoil = false
    transition = false
    similar = false
    forced_transition_flag = false
    ####################### Wake #######################
    for i in 1:length(panels_wake)+1
        if i == 1
            Rs = initialize_wake_first_node(U_nodes[idx_l[end],:], U_nodes[idx_u[end],:], U_nodes[idx_w[1],:], panel_TE, 
                TE_l_isturbulent, TE_u_isturbulent, 
                operatingparameters; init_flag=false)

            Rs_global[idx_w[1],:] .= Rs
        else
            station_residuals = ResidualParamsStation(eltype(U).(sw[i-1:i]), U_nodes[idx_w[i-1],2], U_nodes[idx_w[i-1],1], U_nodes[idx_w[i-1],3], U_nodes[idx_w[i-1]:idx_w[i],4], 
                                                    eltype(U).(operatingparameters.wake_gap[i-1:i]), operatingparameters, upper, 
                                                    turbulent, transition, similar, airfoil, forced_transition_flag)

            Rs_local = station_residuals(U_nodes[idx_w[i],1:3], nothing)
            Rs_global[idx_w[i],:] .= Rs_local
        end
    end

    Rs_global_velocity = get_residual_velocity(U_nodes[:,2], U_nodes[:,4], D, ueinv) #todo: U_nodes seems to have 6 more rows than I want it to have . . .

    if p.operatingparameters.transition_method == 3
        Rs_all = vcat(reshape(Rs_global', :, 1), Rs_global_velocity, Rs_transition)
    else
        Rs_all = vcat(reshape(Rs_global', :, 1), Rs_global_velocity)
    end

    return Rs_all
end

function get_residual_velocity(deltastars, ues, D, ueinv)
    return ues - (ueinv + D*(deltastars .* ues))
end

function create_B_influence_matrix(panels_airfoil, panels_wake, inviscidoutputs)

    numpanels_airfoil = length(panels_airfoil)
    numpanels_wake = length(panels_wake)

    xstar_airfoil = inviscidoutputs.xstar
    ystar_airfoil = inviscidoutputs.ystar
    logr_airfoil = inviscidoutputs.logij
    β_airfoil = inviscidoutputs.β

    B = zeros(numpanels_airfoil+2, numpanels_airfoil+numpanels_wake)
    for i in 1:numpanels_airfoil+1 #airfoil nodes
        if i < numpanels_airfoil+1
            xi = panels_airfoil[i].x1
            yi = panels_airfoil[i].y1
        else
            xi = panels_airfoil[numpanels_airfoil].x2
            yi = panels_airfoil[numpanels_airfoil].y2
        end

        # loop over airfoil panels
        for j in 1:numpanels_airfoil #airfoil panels influencing point i
            # debug = false
            B[i,j] = get_constant_source_streamfunction_influence(xstar_airfoil[i,j], ystar_airfoil[i,j], logr_airfoil[i,j], logr_airfoil[i,j+1], panels_airfoil[j].length, β_airfoil[i,j])#; debug)
        end

        # loop over wake panels
        for k in 1:numpanels_wake
            # define nodes including midpoint
            x1k = panels_wake[k].x1
            x2k = panels_wake[k].x2
            xmk = (x1k + x2k) / 2
            y1k = panels_wake[k].y1
            y2k = panels_wake[k].y2
            ymk = (y1k + y2k) / 2

            if k == numpanels_wake # ghost extension on last point
                x2k = 2 * x2k - xmk
                y2k = 2 * y2k - ymk
                l_l = sqrt((xmk-x1k)^2 + (ymk-y1k)^2)
                l_r = sqrt((x2k-xmk)^2 + (y2k-ymk)^2)
            else
                l_l = l_r = sqrt((x2k-x1k)^2 + (y2k-y1k)^2)/2
            end

            stheta_k_l = (ymk - y1k) / l_l
            ctheta_k_l = (xmk - x1k) / l_l
            stheta_k_r = (y2k - ymk) / l_r
            ctheta_k_r = (x2k - xmk) / l_r

            xstar_ik_l = (xi - x1k) * ctheta_k_l + (yi - y1k) * stheta_k_l
            ystar_ik_l = -(xi - x1k) * stheta_k_l + (yi - y1k) * ctheta_k_l
            xstar_ik_r = (xi - xmk) * ctheta_k_r + (yi - ymk) * stheta_k_r
            ystar_ik_r = -(xi - xmk) * stheta_k_r + (yi - ymk) * ctheta_k_r

            r1k = sqrt((xi-x1k)^2 + (yi-y1k)^2)
            rmk = sqrt((xi-xmk)^2 + (yi-ymk)^2)
            r2k = sqrt((xi-x2k)^2 + (yi-y2k)^2)

            logr1k = iszero(r1k) ? 0.0 : log(r1k)
            logrmk = iszero(rmk) ? 0.0 : log(rmk)
            logr2k = iszero(r2k) ? 0.0 : log(r2k)

            beta_ik_l = get_beta_mfoil(1, 2, xi, yi, x1k, xmk, y1k, ymk) #first two indices just need to be different
            beta_ik_r = get_beta_mfoil(1, 2, xi, yi, xmk, x2k, ymk, y2k)

            # left-panel
            a_l, b_l = get_linear_source_streamfunction_influence(xstar_ik_l, ystar_ik_l, logr1k, logrmk, l_l, beta_ik_l, r1k, rmk)

            # if normal panel, then
            if k > 1
                B[i, numpanels_airfoil+k] += a_l / 2 + b_l
                B[i, numpanels_airfoil+k-1] += a_l / 2
            
            else # if 1st panel
                B[i, numpanels_airfoil+k] += b_l
            end

            # repeat for right panel
            a_r, b_r = get_linear_source_streamfunction_influence(xstar_ik_r, ystar_ik_r, logrmk, logr2k, l_r, beta_ik_r, rmk, r2k)
            
            B[i, numpanels_airfoil+k] += a_r + b_r / 2
            # if normal panel, then
            if k < numpanels_wake
                B[i, numpanels_airfoil+k+1] += b_r / 2
            else # if last panel
                B[i, numpanels_airfoil+k] += b_r / 2
            end
        end
    end

    return B
end

function create_Bprime_influence_matrix(B, inviscidoutputs)

    Bprime = inviscidoutputs.A\B
    Bprime = Bprime[1:end-1,:]

    return Bprime
end

function create_Cgamma_influence_matrix(panels_airfoil, panel_TE, panels_wake)

    blunt_TE = !isnothing(panel_TE)
    numpanels_airfoil = length(panels_airfoil)
    numpanels_wake = length(panels_wake)

    Cgamma = zeros(numpanels_wake+1, numpanels_airfoil+1)
    for k in 1:numpanels_wake+1 # wake nodes
        Vg = zeros(2,numpanels_airfoil+1)

        if k < numpanels_wake+1
            xk = panels_wake[k].x1
            yk = panels_wake[k].y1
        else
            xk = panels_wake[numpanels_wake].x2
            yk = panels_wake[numpanels_wake].y2
        end

        # loop over airfoil panels
        for i in 1:numpanels_airfoil
            x1i = panels_airfoil[i].x1
            x2i = panels_airfoil[i].x2
            y1i = panels_airfoil[i].y1
            y2i = panels_airfoil[i].y2

            li = panels_airfoil[i].length
            stheta_i = panels_airfoil[i].stheta
            ctheta_i = panels_airfoil[i].ctheta

            xstar_ki = (xk - x1i) * ctheta_i + (yk - y1i) * stheta_i
            ystar_ki = -(xk - x1i) * stheta_i + (yk - y1i) * ctheta_i

            r1ki = sqrt((xk-x1i)^2 + (yk-y1i)^2)
            r2ki = sqrt((xk-x2i)^2 + (yk-y2i)^2)

            logr1_ki = iszero(r1ki) ? 0.0 : log(r1ki)
            logr2_ki = iszero(r2ki) ? 0.0 : log(r2ki)

            beta_ki = get_beta_mfoil(1, 2, xk, yk, x1i, x2i, y1i, y2i) #first two indices don't matter

            # get velocity influence from normal airfoil panels 
            a, b = get_linear_vortex_velocity_influence(xstar_ki, ystar_ki, logr1_ki, logr2_ki, beta_ki, li, x1i, y1i, x2i, y2i)
            Vg[:,i] .+= a
            Vg[:,i+1] .+= b
        end

        # now TE panel
        if blunt_TE
            x1_TE = panel_TE.x1
            y1_TE = panel_TE.y1
            x2_TE = panel_TE.x2
            y2_TE = panel_TE.y2

            l_TE = panel_TE.length
            stheta_TE = panel_TE.stheta
            ctheta_TE = panel_TE.ctheta

            xstar_kTE = (xk - x1_TE) * ctheta_TE + (yk - y1_TE) * stheta_TE
            ystar_kTE = -(xk - x1_TE) * stheta_TE + (yk - y1_TE) * ctheta_TE

            r1_kTE = sqrt((xk - x1_TE)^2 + (yk - y1_TE)^2)
            r2_kTE = sqrt((xk - x2_TE)^2 + (yk - y2_TE)^2)

            logr1_kTE = iszero(r1_kTE) ? 0.0 : log(r1_kTE)
            logr2_kTE = iszero(r2_kTE) ? 0.0 : log(r2_kTE)

            beta_kTE = get_beta_mfoil(1,2, xk, yk, x1_TE, x2_TE, y1_TE, y2_TE) #indices are meaningless as long as they are different

            # source
            a = get_constant_source_velocity_influence(xstar_kTE, ystar_kTE, l_TE, logr1_kTE, logr2_kTE, beta_kTE, panel_TE.x1, panel_TE.y1, panel_TE.x2, panel_TE.y2)
            f = a * panel_TE.tcp / 2
            Vg[:,1] .-= f
            Vg[:,numpanels_airfoil+1] .+= f

            #vortex
            a, b = get_linear_vortex_velocity_influence(xstar_kTE, ystar_kTE, logr1_kTE, logr2_kTE, beta_kTE, panel_TE.length, panel_TE.x1, panel_TE.y1, panel_TE.x2, panel_TE.y2)
            g = (a + b) * panel_TE.tdp / 2
            Vg[:,1] .+= g
            Vg[:,numpanels_airfoil+1] .-= g
        end

        if k < numpanels_wake+1
            Cgamma[k,:] = Vg[1,:] * panels_wake[k].tdir1[1] + Vg[2,:] * panels_wake[k].tdir1[2]
        else
            Cgamma[k,:] = Vg[1,:] * panels_wake[k-1].tdir2[1] + Vg[2,:] * panels_wake[k-1].tdir2[2]
        end
    end

    return Cgamma
end

function create_Csigma_influence_matrix(panels_airfoil, panels_wake)

    numpanels_airfoil = length(panels_airfoil)
    numpanels_wake = length(panels_wake)

    Csigma = zeros(numpanels_wake+1, numpanels_airfoil+numpanels_wake)
    #loop over wake nodes (control points)
    for k in 1:numpanels_wake+1

        if k < numpanels_wake+1
            xk = panels_wake[k].x1
            yk = panels_wake[k].y1
            wake_tdir = panels_wake[k].tdir1
        else
            xk = panels_wake[numpanels_wake].x2
            yk = panels_wake[numpanels_wake].y2
            wake_tdir = panels_wake[k-1].tdir2
        end
        
        # loop over airfoil source panels
        for i in 1:numpanels_airfoil
            if !(i==1 && k==1) && !(i==numpanels_airfoil && k==1) #first wake node is affected differently by first and last panels on airfoil #todo: include what that is
                x1i = panels_airfoil[i].x1
                x2i = panels_airfoil[i].x2
                y1i = panels_airfoil[i].y1
                y2i = panels_airfoil[i].y2

                li = panels_airfoil[i].length
                stheta_i = panels_airfoil[i].stheta
                ctheta_i = panels_airfoil[i].ctheta

                xstar_ki = (xk - x1i) * ctheta_i + (yk - y1i) * stheta_i
                ystar_ki = -(xk - x1i) * stheta_i + (yk - y1i) * ctheta_i

                r1ki = sqrt((xk-x1i)^2 + (yk-y1i)^2)
                r2ki = sqrt((xk-x2i)^2 + (yk-y2i)^2)

                logr1_ki = iszero(r1ki) ? 0.0 : log(r1ki)
                logr2_ki = iszero(r2ki) ? 0.0 : log(r2ki)

                beta_ki = get_beta_mfoil(k, i, xk, yk, x1i, x2i, y1i, y2i)

                Csigma[k,i] = get_constant_source_velocity_influence(xstar_ki, ystar_ki, li, logr1_ki, logr2_ki, beta_ki, x1i, y1i, x2i, y2i; vdir=wake_tdir) 
            end
        end 

        # loop over linear source wake panels
        for j in 1:numpanels_wake+1
            # define nodes including midpoint
            if j == 1
                x1j = panels_wake[j].x1
                xmj = panels_wake[j].x1
                x2j = panels_wake[j].xbar
                y1j = panels_wake[j].y1
                ymj = panels_wake[j].y1
                y2j = panels_wake[j].ybar
            elseif j == numpanels_wake+1
                x1j = panels_wake[j-1].xbar
                xmj = panels_wake[j-1].x2
                x2j = 2 * xmj - x1j
                y1j = panels_wake[j-1].ybar
                ymj = panels_wake[j-1].y2
                y2j = 2 * ymj - y1j
            else
                x1j = panels_wake[j-1].xbar
                xmj = panels_wake[j].x1
                x2j = panels_wake[j].xbar
                y1j = panels_wake[j-1].ybar
                ymj = panels_wake[j].y1
                y2j = panels_wake[j].ybar
            end

            l_l = sqrt((xmj-x1j)^2 + (ymj-y1j)^2)
            l_r = sqrt((x2j-xmj)^2 + (y2j-ymj)^2)

            stheta_j_l = (ymj - y1j) / l_l
            ctheta_j_l = (xmj - x1j) / l_l
            stheta_j_r = (y2j - ymj) / l_r
            ctheta_j_r = (x2j - xmj) / l_r

            xstar_kj_l = (xk - x1j) * ctheta_j_l + (yk - y1j) * stheta_j_l
            ystar_kj_l = -(xk - x1j) * stheta_j_l + (yk - y1j) * ctheta_j_l
            xstar_kj_r = (xk - xmj) * ctheta_j_r + (yk - ymj) * stheta_j_r
            ystar_kj_r = -(xk - xmj) * stheta_j_r + (yk - ymj) * ctheta_j_r

            r1kj = sqrt((xk-x1j)^2 + (yk-y1j)^2)
            rmkj = sqrt((xk-xmj)^2 + (yk-ymj)^2)
            r2kj = sqrt((xk-x2j)^2 + (yk-y2j)^2)

            logr1kj = iszero(r1kj) ? 0.0 : log(r1kj)
            logrmkj = iszero(rmkj) ? 0.0 : log(rmkj)
            logr2kj = iszero(r2kj) ? 0.0 : log(r2kj)

            beta_kj_l = get_beta_mfoil(1, 2, xk, yk, x1j, xmj, y1j, ymj) #first two indices just need to be different
            beta_kj_r = get_beta_mfoil(1, 2, xk, yk, xmj, x2j, ymj, y2j)

            if k == j #control point on same panel
                if j==1
                    l_lower = panels_airfoil[1].length
                    l_upper = panels_airfoil[end].length

                    Csigma[k,1] += 1/(2pi) * (log(l_lower/l_r) + 1)
                    Csigma[k,numpanels_airfoil] += 1/(2pi) * (log(l_upper/l_r) + 1)
                    Csigma[k,numpanels_airfoil+1] -= 1/(2pi)
                elseif j < numpanels_wake+1
                    aa = 1 / (4pi) * log(l_l/l_r)

                    Csigma[k,numpanels_airfoil+j-1] += aa + 1/(2pi)
                    Csigma[k,numpanels_airfoil+j] += aa - 1/(2pi)
                end
            else
                if j==1 #only right-panel
                    a, b = get_linear_source_velocity_influence(xstar_kj_r, ystar_kj_r, logrmkj, logr2kj, beta_kj_r, l_r, xmj, ymj, x2j, y2j; vdir=wake_tdir)

                    Csigma[k,numpanels_airfoil+1] += b
                    Csigma[k,1] += a
                    Csigma[k,numpanels_airfoil] += a
                elseif j==numpanels_wake+1 #last point has a constant source ghost extension
                    l_end = sqrt((x2j-x1j)^2 + (y2j-y1j)^2)
                    stheta_j_end = (y2j - y1j) / l_end
                    ctheta_j_end = (x2j - x1j) / l_end
                    xstar_kj_end = (xk - x1j) * ctheta_j_end + (yk - y1j) * stheta_j_end
                    ystar_kj_end = -(xk - x1j) * stheta_j_end + (yk - y1j) * ctheta_j_end
                    beta_kj_end = get_beta_mfoil(1, 2, xk, yk, x1j, x2j, y1j, y2j)
                    
                    a = get_constant_source_velocity_influence(xstar_kj_end, ystar_kj_end, l_end, logr1kj, logr2kj, beta_kj_end, x1j, y1j, x2j, y2j; vdir=wake_tdir)

                    Csigma[k,numpanels_airfoil+numpanels_wake] += a
                else
                    a1, b1 = get_linear_source_velocity_influence(xstar_kj_l, ystar_kj_l, logr1kj, logrmkj, beta_kj_l, l_l, x1j, y1j, xmj, ymj; vdir=wake_tdir) # left half-panel
                    a2, b2 = get_linear_source_velocity_influence(xstar_kj_r, ystar_kj_r, logrmkj, logr2kj, beta_kj_r, l_r, xmj, ymj, x2j, y2j; vdir=wake_tdir) # right half-panel

                    Csigma[k,numpanels_airfoil+j-1] += a1 + b1/2
                    Csigma[k,numpanels_airfoil+j] += a2/2 + b2
                end
            end
        end
    end

    return Csigma
end

function create_Dwake_influence_matrix(Bprime, Cgamma, Csigma)

    Dwake = Cgamma * Bprime + Csigma
    Dwake[1,:] = Bprime[end,:]

    return Dwake
end

function create_Dprime_influence_matrix(panels_airfoil, panels_wake, d_airfoil)
    
    numpanels_airfoil = length(panels_airfoil)
    numpanels_wake = length(panels_wake)

    # get s
    s_airfoil, s_wake, _ = get_s_length_vectors(panels_airfoil, panels_wake)

    return create_Dprime_influence_matrix(s_airfoil, s_wake, d_airfoil, nothing)
end

function create_Dprime_influence_matrix(s_airfoil, s_wake, d_airfoil, d)

    numpanels_airfoil = length(s_airfoil) - 1
    numpanels_wake = length(s_wake) - 1

    Dprime = spzeros(numpanels_airfoil+numpanels_wake, numpanels_airfoil+numpanels_wake+2)
    #loop through airfoil directions
    for i in 1:numpanels_airfoil
        ds = s_airfoil[i+1] - s_airfoil[i]
        Dprime[i,i] = -d_airfoil[i]/ds
        Dprime[i,i+1] = d_airfoil[i+1]/ds
    end

    #loop through wake panels
    for k in 1:numpanels_wake
        ds = s_wake[k+1] - s_wake[k]
        Dprime[numpanels_airfoil+k,numpanels_airfoil + 1 + k] = -1/ds #direction of wake will always be positive
        Dprime[numpanels_airfoil+k,numpanels_airfoil + 1 + k+1] = 1/ds
    end

    return Dprime
end

function get_d_directional_vectors(numpanels_wake, inviscidoutputs; alpha_idx=1)

    d_airfoil = vcat(sign.(inviscidoutputs.solution[:, alpha_idx]))
    d_wake = ones(numpanels_wake+1)
    d = vcat(d_airfoil, d_wake) 

    return d_airfoil, d_wake, d
end

function get_s_length_vectors(panels_airfoil, panels_wake, s_stag=0.0)

    s_airfoil = vcat(0.0, cumsum([p.length for p in panels_airfoil]))
    s_wake = vcat(s_airfoil[end], s_airfoil[end] .+ cumsum([p.length for p in panels_wake]))

    return s_airfoil, s_wake .- s_stag, s_wake
end

function create_D_influence_matrix(panels_airfoil, panel_TE, panels_wake, inviscidoutputs; alpha_idx=1)

    B = create_B_influence_matrix(panels_airfoil, panels_wake, inviscidoutputs)

    Bprime = create_Bprime_influence_matrix(B, inviscidoutputs)

    Cgamma = create_Cgamma_influence_matrix(panels_airfoil, panel_TE, panels_wake)

    Csigma = create_Csigma_influence_matrix(panels_airfoil, panels_wake)

    Dwake = create_Dwake_influence_matrix(Bprime, Cgamma, Csigma)

    d_airfoil, _, d = get_d_directional_vectors(length(panels_wake), inviscidoutputs; alpha_idx)

    Dprime = create_Dprime_influence_matrix(panels_airfoil, panels_wake, d_airfoil)

    d_diag = spdiagm(d)

    D = d_diag * vcat(Bprime, Dwake) * Dprime

    return D, d, Bprime, Dwake
end

function recalculate_D(Bprime, Dwake, d_airfoil, d, s_airfoil, s_wake)

    Dprime = create_Dprime_influence_matrix(s_airfoil, s_wake, d_airfoil, d)
    # Dprime = create_Dprime_influence_matrix(panels_airfoil, panels_wake, d_airfoil)
    d_diag = spdiagm(d)

    return d_diag * vcat(Bprime, Dwake) * Dprime
end

function build_wake(operatingparameters, panels_airfoil, panel_TE, s, γref, sinalpha, cosalpha)

    numpanels_wake = operatingparameters.numpanels_wake
    wake_nodes = zeros(2, numpanels_wake+1)
    wake_tangential_directions = zeros(2, numpanels_wake+1)

    # build out wake points based on streamline
    ds1 = max((s[2] - s[1] + s[end] - s[end-1]) / 2, 1e-3) #first nominal wake panel size
    wake_length = operatingparameters.wake_length * operatingparameters.chord
    sv = geometric_spaced_points(ds1, wake_length, numpanels_wake+1)
    # sv = geometric_points(ds1, wake_length, numpanels_wake+1)

    # first point just behind the TE
    xy1 = [panels_airfoil[1].x1, panels_airfoil[1].y1]
    xyN = [panels_airfoil[end].x2, panels_airfoil[end].y2]
    n_TE = xyN - xy1 # normal vector
    t_TE = [n_TE[2], -n_TE[1]] #tangential vector
    TE_midpoint = (xy1 + xyN) / 2
    wake_nodes[:,1] = TE_midpoint + 1e-5 * t_TE * operatingparameters.chord

    γs = γref * [cosalpha, sinalpha]

    Vinf = operatingparameters.Vinf

    # loop through rest of the wake nodes
    for i in 1:numpanels_wake
        # velocity at first node
        v1 = get_invisicd_velocity_from_gammas(panels_airfoil, γs, wake_nodes[:,i], panel_TE, operatingparameters; cosalpha, sinalpha)
        v1 = v1 / norm(v1) #normalize
        wake_tangential_directions[:,i] = v1 #save for later
        wake_nodes[:,i+1] = wake_nodes[:,i] + (sv[i+1] - sv[i]) * v1 # predict next node location

        # velocity at second node - corrector step
        v2 = get_invisicd_velocity_from_gammas(panels_airfoil, γs, wake_nodes[:,i+1], panel_TE, operatingparameters; cosalpha, sinalpha)

        v2 = v2 / norm(v2) 
        wake_tangential_directions[:,i+1] = v2 
        wake_nodes[:,i+1] = wake_nodes[:,i] + (sv[i+1] - sv[i]) * (v1 + v2) / 2 # actual step is an average of v1 and v2
    end

    ue_wake_ref = zeros(numpanels_wake+1, 2)
    # once wake nodes are defined, re-find the velocity at each node
    for i in 1:numpanels_wake+1
        v_0 = get_invisicd_velocity_from_gammas(panels_airfoil, γref[:,1], wake_nodes[:,i], panel_TE, operatingparameters; cosalpha=1.0, sinalpha=0.0)
        v_90 = get_invisicd_velocity_from_gammas(panels_airfoil, γref[:,2], wake_nodes[:,i], panel_TE, operatingparameters; cosalpha=0.0, sinalpha=1.0)

        ue_wake_ref[i,1] = v_0' * wake_tangential_directions[:,i]
        ue_wake_ref[i,2] = v_90' * wake_tangential_directions[:,i]
    end

    # fill in the rest of the stuff I need for wake panels
    xwake = wake_nodes[1,:]
    ywake = wake_nodes[2,:]

    xbars = (xwake[1:end-1] .+ xwake[2:end])/2
    ybars = (ywake[1:end-1] .+ ywake[2:end])/2
    lengths = @. sqrt((xwake[2:end] - xwake[1:end-1]) ^ 2 + (ywake[2:end] - ywake[1:end-1]) ^ 2)

    sthetas = @. (ywake[2:end] - ywake[1:end-1]) / lengths
    cthetas = @. (xwake[2:end] - xwake[1:end-1]) / lengths

    tdir = [wake_tangential_directions[:,i] for i in 1:size(wake_tangential_directions, 2)]

    panels_wake = WakePanel{eltype(xwake)}.(xwake[1:end-1], ywake[1:end-1], xwake[2:end], ywake[2:end], xbars, ybars, lengths, sthetas, cthetas, tdir[1:end-1], tdir[2:end]) 

    wake_gap = zeros(numpanels_wake+1)
    wake_node_distance = cumsum(lengths)
    f_length = 2.5
    dtdx = min(max(panel_TE.dtdx, -3/f_length), 3/f_length) #clip TE thickness slope
    Lw = f_length * panel_TE.hTE
    wake_gap[1] = panel_TE.hTE
    for i in 2:numpanels_wake+1
        xib = wake_node_distance[i-1] / Lw
        if xib < 1
            wake_gap[i] = panel_TE.hTE * (1 + (2 + f_length * dtdx) * xib) * (1 - xib)^2
        end
    end

    operatingparameters.wake_gap[:] .= wake_gap

    return panels_wake, ue_wake_ref
end

function stagnation_point_move!(p, ue) #ue will be for all nodes, but should be all positive, as if it's the boundary layer perspective

    (; s, idx_u, idx_l, idx_w, operatingparameters, 
        D, ueinv, idx_transition_u, idx_transition_l, 
        panels_airfoil, panel_TE, panels_wake, inviscidoutputs, d, Bprime, Dwake, swake0) = p

    # for specifics, see MFOIL 2243, and find a specific example of it happening to really get it
    newpanel_flag = false

    if (ue[idx_u[1]] > 0 && any(ue[idx_u] .< 0)) || (ue[idx_l[1]] > 0 && any(ue[idx_l] .< 0))
        # negative velocity, but not the first one, so the solver is updating a downstream velocity state incorrectly (if it should be negative, so should the preceding states)
        error("One or more surfaces contains ")
    end

    if ue[idx_u[1]] < 0 # move stagnation point up

        panel_steps_up = findfirst(x->x > 0, ue[idx_u[1]:idx_u[end]])-1
        newpanel_idx = idx_u[1]+panel_steps_up
        for i in idx_u[1]:newpanel_idx-1 # often this will just be looping over one panel - but just in case, it could be more
            ue[i] = -ue[i]
            d[i] = -d[i]
            p.ueinv[i] *= -1
        end

        # shift idx vectors
        idx_u[1] += panel_steps_up
        idx_l[1] += panel_steps_up

        #shift transition idx
        idx_transition_u[1] -= panel_steps_up
        idx_transition_l[1] += panel_steps_up

        newpanel_flag = true

    elseif ue[idx_l[1]] < 0 # move stagnation panel down

        panel_steps_down = findfirst(x->x > 0, ue[idx_l[1]:-1:idx_l[end]])-1
        newpanel_idx = idx_l[1]-panel_steps_down
        for i in idx_l[1]:-1:newpanel_idx+1 # often this will just be looping over one panel - but just in case, it could be more
            ue[i] = -ue[i]
            d[i] = -d[i]
            p.ueinv[i] *= -1
        end

        # shift idx vectors
        idx_u[1] -= panel_steps_down
        idx_l[1] -= panel_steps_down

        #shift transition idx
        idx_transition_u[1] += panel_steps_down
        idx_transition_l[1] -= panel_steps_down
        
        newpanel_flag = true
    end # else panel is the same

    # ue[findall(x->x<0, ue)] .*= -1 #shouldn't need this anymore, right?

    # recalculate exact stagnation point
    su, sl, sw, s_stag = split_s(s, ue, idx_u, idx_l, idx_w, swake0)

    if newpanel_flag # recalculate D
        d_airfoil = d[1:idx_u[end]]
        Dnew = recalculate_D(Bprime, Dwake, d_airfoil, d, s, sw)
        D[:] = Dnew[:]
    end

    return su, sl, sw, s_stag, ue #convenient to not solve again in parent function
end

function resolve_transition_location!(p, U; upper=true)

    (; operatingparameters, D, ueinv, d, turbulent_nodes) = p

    if !upper
        idx_surface = p.idx_l
        idx_nodes = idx_surface[1]:-1:idx_surface[end]
        # s = p.sl
        idx_transition = p.idx_transition_l
    else
        idx_surface = p.idx_u
        idx_nodes = idx_surface[1]:idx_surface[end]
        # s = p.su
        idx_transition = p.idx_transition_u
    end

    su, sl, _, _ = split_s(p.s, U[:,4], p.idx_u, p.idx_l, p.idx_w, p.swake0)

    if !upper
        s = sl
    else
        s = su
    end

    amp_states = U[idx_nodes, 3] # save for later
    
    ### resolve 3rd state along boundary layer
    idx_turb = march_amplification!(s, U, turbulent_nodes, idx_nodes, operatingparameters; upper)

    if idx_turb == idx_transition[1]
        U[idx_nodes,3] .= amp_states #don't accept changes if don't need to
        return
    end

    if operatingparameters.verbose
        surface = upper ? "upper" : "lower"
        println("Update ", surface, " transition: last lam [$(idx_transition[1]-1)]->[$(idx_turb-1)]")
    end

    if idx_turb < idx_transition[1] # transition earlier than before --> fill in new turbulent states

        # get ct0 turbulent initial value
        ct0 = get_ct_transition(U[idx_nodes[idx_turb-1],:], operatingparameters; turbulent=true)

        if idx_transition[1] < length(s)
            ct1 = U[idx_nodes[idx_transition[1]],3]
        else
            ct1 = ct0
        end

        dx = s[min(idx_transition[1], length(s))] - s[idx_turb]

        for i in (idx_turb):min(idx_transition[1]-1, length(idx_nodes))
            if dx == 0 || i == (idx_turb)
                f = 0
            elseif idx_turb == idx_transition[1]-1
                f = 1
            else
                f = (s[i] - s[idx_turb]) / dx
            end

            U[idx_nodes[i],3] = ct0 + f * (ct1 - ct0)
            turbulent_nodes[idx_nodes[i]] = true
        end

    else # transition is later ---> leave alone. Just change the turbulent nodes vector
        for i in idx_transition[1]:idx_turb-1
            turbulent_nodes[idx_nodes[i]] = false
        end
    end

    # move xt states to mid-panel on the new panel
    if p.operatingparameters.transition_method == 3
        if idx_turb == length(s)+1
            if upper 
                U[end-2] = s[idx_turb-1] #upper xt
                U[end-1] = U[idx_nodes[idx_turb-1], 1] #upper theta_t
                U[end] = U[idx_nodes[idx_turb-1], 2] #upper delta_t
            else
                U[end-5] = s[idx_turb-1] #lower xt
                U[end-4] = U[idx_nodes[idx_turb-1], 1] #upper theta_t
                U[end-3] = U[idx_nodes[idx_turb-1], 2] #upper delta_t
            end
        else
            s1 = s[idx_turb-1]
            s2 = s[idx_turb]
            if upper 
                U[end-2] = (s1+s2)/2 #upper xt
                U[end-1] = (U[idx_nodes[idx_turb-1],1] + U[idx_nodes[idx_turb],1])/2 #upper theta_t
                U[end] = (U[idx_nodes[idx_turb-1],2] + U[idx_nodes[idx_turb],2])/2 #upper delta_t
            else
                U[end-5] = (s1+s2)/2 #lower xt
                U[end-4] = (U[idx_nodes[idx_turb-1],1] + U[idx_nodes[idx_turb],1])/2 #upper theta_t
                U[end-3] = (U[idx_nodes[idx_turb-1],2] + U[idx_nodes[idx_turb],2])/2 #upper delta_t
            end
        end
    end

    # update transition location
    idx_transition[1] = idx_turb
end

function march_amplification!(s, U, turbulent_nodes, idx_nodes, operatingparameters; upper=true)

    turbulent = false

    U[idx_nodes[1],3] = 0.0 # start at zero

    sft = upper ? operatingparameters.sft_u[1] : operatingparameters.sft_l[1]

    i = 2
    while i <= length(idx_nodes)

        resolve_transition_params = ResolveTransitionParams([s[i-1], s[i]], U[idx_nodes[i-1],:], U[idx_nodes[i],:], operatingparameters, turbulent)

        initial_values = turbulent_nodes[idx_nodes[i]] ? [U[idx_nodes[i-1],3]*1.01] : [U[idx_nodes[i],3]]

        debug = false

        solution, converged = local_newton_amp(resolve_transition_params, initial_values; tol=operatingparameters.etol, max_iter=operatingparameters.maxiters, debug)

        if !converged
            solution = initial_values
        end

        n2 = solution[1]

        #todo: eventually implement the forced transition
        # check for transition
        if n2 > operatingparameters.ncrit || s[i] > sft
            idx_turb = i
            return idx_turb
        else # Assign solution to U2
            U[idx_nodes[i], 3] = n2
        end

        i += 1
    end

    idx_turb = length(idx_nodes) + 1 # if never transitions

    return idx_turb
end

struct ResolveTransitionParams{TF}
    xs::AbstractVector{TF}
    U1::AbstractArray{TF}
    U2::AbstractArray{TF}
    params::OperatingParameters{TF}
    turbulent::Bool
end

function (p::ResolveTransitionParams)(states, params)

    U2 = eltype(states).(p.U2)
    U2[3] = states[1]

    R_amp = residual_amplification(p.xs, p.U1, U2, p.params)

    return R_amp
end

(p::ResolveTransitionParams)(states) = (p::ResolveTransitionParams)(states, nothing)
