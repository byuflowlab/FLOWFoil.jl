
############ Combined problem ##############
function initialize_boundary_layer(s, operatingparameters, ues=nothing, U_stag_wake=nothing; 
                                    upper=true, airfoil=true, plot_bool=false, kwargs...)

    (;rho0, mu0, ncrit) = operatingparameters

    npanels = length(s) - 1

    # transition definitions
    idx_transition = npanels+2 # default no transition
    forced_transition_flag = false
    if airfoil
        # forced_transition_flag = upper ? operatingparameters.forced_u : operatingparameters.forced_l
        forced_transition_flag = false #with mfoil's approach, assume false at first

        if upper
            sft = operatingparameters.sft_u[1]
        else
            sft = operatingparameters.sft_l[1]
        end

        #determine idx_transition if forced
        # if forced_transition_flag
        #     idx_transition = findfirst(x->x>sft, s)
        # end
    else
        forced_transition_flag = false #no transition in the wake
    end

    # if !forced_transition_flag
        deltas = zeros(npanels+1)
        thetas = zeros(npanels+1) 
        state3s = zeros(npanels+1) 
        if isnothing(ues)
            ues = ones(npanels+1)
        end
    # else
        # deltas = zeros(npanels+2)
        # thetas = zeros(npanels+2) 
        # state3s = zeros(npanels+2) 
        # if isnothing(ues) #? do I actually ever do this?
        #     ues = ones(npanels+2)
        # else
            # redefine ues vector to add in transition state I haven't determined yet, linearly interpolating for now
                # I'm ok redefining this to a new length vector since I only do it this once for each airfoil surface
            # ue_t = FLOWMath.linear(s[idx_transition-1:idx_transition], ues[idx_transition-1:idx_transition], sft)
            # ues = vcat(ues[1:idx_transition-1], ue_t, ues[idx_transition:end])
        # end
    # end

    transition = false
    theta_t = 0.0
    delta_t = 0.0

    st = forced_transition_flag ? sft : s[end]

    if airfoil
        if s[1] > 1e-8
            K = ues[1] / s[1]
        else
            K = ues[2] / s[2]
        end
        θ0 = sqrt(0.45 * mu0 / (6 * rho0 * K))
        δstar0 = 2.2 * θ0
        ñ0 = 0.0

        turbulent = false #start laminar in airfoil case

        # solve first node from stagnation state
        if s[1] < 1e-8 * s[end]
            hitstag = true
        else
            hitstag = false
        end
        
        sstag = 1e-6
        if !hitstag
            ue_stag = ues[1]/s[1] * sstag
        else
            ue_stag = ues[2]/s[2] * sstag
        end

        U_stag = station_stagnation(δstar0, θ0, ñ0, [ue_stag, ue_stag], [sstag, sstag], operatingparameters; upper)

        if !hitstag
            thetas[1] = U_stag[1]
            deltas[1] = U_stag[2]
            state3s[1] = U_stag[3]
    
            i = 1
        else
            thetas[1] = U_stag[1]
            deltas[1] = U_stag[2]
            state3s[1] = U_stag[3]
            
            thetas[2] = U_stag[1]
            deltas[2] = U_stag[2]
            state3s[2] = U_stag[3]
    
            i = 2
        end

        if sft < s[2]
            turbulent = true # transition triggered at the beginning of the boundary layer, can proceed with turbulent flow
        end

        wake_gap = zeros(npanels+1)
    else #wake 
        @assert !isnothing(U_stag_wake) "For the wake boundary layer, the U_stag_wake argument must be given."
        U_stag = U_stag_wake
        turbulent = true # start turbulent in the wake
        thetas[1] = U_stag[1]
        deltas[1] = U_stag[2]
        state3s[1] = U_stag[3]

        i = 1
        wake_gap = operatingparameters.wake_gap
    end

    debug = false
    while i <= npanels # use a while loop to be able to navigate it easier
        
        if operatingparameters.transition_method == 1 || (operatingparameters.transition_method == 2 && !forced_transition_flag)
            
            if airfoil && !turbulent && s[i+1] > sft
                transition = true
                turbulent = true
                forced_transition_flag = true
                idx_transition = i+1
            end
            
            if !turbulent
                δstar_temp, θ_temp, state3_temp, ue_temp = station_laminar(deltas[i], thetas[i], state3s[i], ues[i:i+1], s[i:i+1], operatingparameters; upper, debug)
            else
                if transition
                    δstar_temp, θ_temp, state3_temp, ue_temp, st, theta_t, delta_t = station_transition(deltas[i], thetas[i], state3s[i], ues[i:i+1], s[i:i+1], operatingparameters; upper, debug, forced_transition_flag)
                    transition = false
                    forced_transition_flag = false
                else
                    δstar_temp, θ_temp, state3_temp, ue_temp = station_turbulent(deltas[i], thetas[i], state3s[i], ues[i:i+1], s[i:i+1], operatingparameters, wake_gap[i:i+1]; upper, airfoil, debug)
                end
            end

            if !turbulent && state3_temp > ncrit #free transition
                transition = true
                turbulent = true
                idx_transition = i+1
                i -= 1 # repeat this node (it should be unchanged by the end of the while loop)
            else
                deltas[i+1] = δstar_temp
                thetas[i+1] = θ_temp
                state3s[i+1] = state3_temp
                ues[i+1] = ue_temp
            end
        elseif (operatingparameters.transition_method==2 && forced_transition_flag) #forced transition #todo: in geometry when I define sft define forced or not based on xft
            # @infiltrate
            if i < idx_transition #laminar flow
                turbulent = false
                δstar_temp, θ_temp, state3_temp, ue_temp = station_laminar(deltas[i], thetas[i], state3s[i], ues[i:i+1], s[i:i+1], operatingparameters; upper, debug)
            else #turbulent flow
                turbulent = true
                δstar_temp, θ_temp, state3_temp, ue_temp = station_turbulent(deltas[i], thetas[i], state3s[i], ues[i:i+1], s[i:i+1], operatingparameters, wake_gap[i:i+1]; upper, airfoil, debug)
            end

            deltas[i+1] = δstar_temp
            thetas[i+1] = θ_temp
            state3s[i+1] = state3_temp
            ues[i+1] = ue_temp
        end
        i += 1
    end

    return deltas, thetas, state3s, ues, idx_transition, st, theta_t, delta_t
end

struct ResidualParamsStation{TF}
    xs::AbstractVector{TF}
    δstar0::TF
    θ0::TF
    state30::TF
    ues::Vector{TF}
    wake_gap::AbstractVector{TF}    # wake TE gap to be removed from delta* in the wake residuals
    params::OperatingParameters
    upper::Bool                     # on upper or lower surface
    turbulent::Bool                 # turbulent or laminar flow
    transition::Bool                # transition station
    similar::Bool                   # stagnation station
    airfoil::Bool                   # airfoil or wake
    forced_transition_flag::Bool    # in case of transition, is it forced or not?
end

function (p::ResidualParamsStation)(states, params)

    if length(states) == 3
        if p.similar
            deltastars = [states[2], states[2]]
            thetas = [states[1], states[1]]
            state3s = [states[3], states[3]]
        else
            deltastars = [p.δstar0, states[2]]
            thetas = [p.θ0, states[1]]
            state3s = [p.state30, states[3]]
        end
        ues = p.ues
    elseif length(states) == 4 # inverse method needs ue as a state
        if p.similar
            deltastars = [states[2], states[2]]
            thetas = [states[1], states[1]]
            state3s = [states[3], states[3]]
            ues = [states[4], states[4]]
        else
            deltastars = [p.δstar0, states[2]]
            thetas = [p.θ0, states[1]]
            state3s = [p.state30, states[3]]
            ues = [p.ues[1], states[4]]
        end
    elseif length(states) == 6 && p.params.transition_method == 3
            deltastars = [p.δstar0, states[2], states[6]]
            thetas = [p.θ0, states[1], states[5]]
            state3s = [p.state30, states[3]]
            ues = p.ues
            xs = vcat(p.xs, states[4])

    elseif length(states) == 7 && p.params.transition_method == 3
            deltastars = [p.δstar0, states[2], states[7]]
            thetas = [p.θ0, states[1], states[6]]
            state3s = [p.state30, states[3]]
            ues = [p.ues[1], states[4]]
            xs = vcat(p.xs, states[5])
    else
        error("Error: states needs to be of length 3 or 4")
    end

    debug = false
    if p.transition
        if p.params.transition_method == 1
            R = residuals_transition(p.xs, deltastars, thetas, state3s, ues, p.params; upper=p.upper, debug, forced_transition_flag=p.forced_transition_flag)
        # elseif p.params.transition_method == 2
        #     R = residuals_transition_test2(p.xs, deltastars, thetas, state3s, ues, p.params; upper=p.upper, debug)
        # elseif p.params.transition_method == 3
        #     R = residuals_transition_test3(xs, deltastars, thetas, state3s, ues, p.params; upper=p.upper, debug)
        end
    else
        R = residuals_together_local(p.xs, deltastars, thetas, state3s, ues, p.params, p.wake_gap; airfoil=p.airfoil, turbulent=p.turbulent, upper=p.upper, similar=p.similar)
    end
            # @infiltrate
    return R
end

(p::ResidualParamsStation)(states) = (p::ResidualParamsStation)(states, nothing)

struct IsOutOfDomain{TF}
    lower_bounds::AbstractVector{TF}
end

function (iood::IsOutOfDomain)(states, params)
    return any(states .< iood.lower_bounds)
end

function station_turbulent(deltastar0, theta0, ct120, ues, xs, params, wake_gap=zeros(2); plot_bool=true, upper=true, airfoil=true, debug=false)

    (; etol, verbose) = params
    
    initial_values = vcat(theta0, deltastar0, ct120, ues[2])

    turbulent = true
    transition = false
    similar = false
    forced_transition_flag = transition

    params_solver = ResidualParamsStation(xs, deltastar0, theta0, ct120, ues, wake_gap, params, upper, turbulent, transition, similar, airfoil, forced_transition_flag)

    if params.use_native_solver_local
        # solution, converged = local_levenberg_marquardt(params_solver, initial_values; tol=etol, max_iter=params.maxiters)
        solution, converged = local_newton(params_solver, initial_values; tol=etol, max_iter=params.maxiters, debug)

        # if unconverged, extrapolate instead
        if !converged
            solution = initial_values
            if airfoil
                solution[1] = theta0 * sqrt(xs[2]/xs[1])
                solution[2] = deltastar0 * sqrt(xs[2]/xs[1])
            else
                rlen = (xs[2]-xs[1]) / (10 * theta0)
                solution[1] = (theta0 + rlen * deltastar0) / (1 + rlen)
            end
            if params.verbose; @warn "BL initialization station not converged. Extrapolating values..."; end
        end
    else

        prob = NonlinearSolve.NonlinearProblem(params_solver, initial_values)
        isoutofdomain = IsOutOfDomain(fill(1e-6, 3)) #todo: generalize better or make a user input?
        sol = NonlinearSolve.solve(prob, NonlinearSolve.LevenbergMarquardt(), 
                                    abstol=etol, 
                                    maxiters=params.maxiters,
                                    isoutofdomain=isoutofdomain)
                                    
        if sol.retcode != NonlinearSolve.ReturnCode.Success
            sol2 = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(), 
            abstol=etol, 
            maxiters=params.maxiters,
            isoutofdomain=isoutofdomain)
            if sol2.retcode == NonlinearSolve.ReturnCode.Success
                sol = sol2
            end
        end
        solution = sol.u #todo: extrapolate if not converging
    end

    theta_2 = solution[1]
    delta_2 = solution[2]
    ct_2 = solution[3]
    ue_2 = solution[4]

    return delta_2, theta_2, ct_2, ue_2
end

function station_laminar(deltastar0, theta0, n0, ues, xs, params; plot_bool=false, upper=true, debug=false)

    (; etol, verbose) = params
    
    initial_values = vcat(theta0, deltastar0, n0, ues[2])

    turbulent = false
    transition = false
    similar = false
    airfoil = true #for laminar, will always be true
    forced_transition_flag = transition

    wake_gap = zeros(2)
    params_solver = ResidualParamsStation(xs, deltastar0, theta0, n0, ues, wake_gap, params, upper, turbulent, transition, similar, airfoil, forced_transition_flag)

    if params.use_native_solver_local
        # solution, converged = local_levenberg_marquardt(params_solver, initial_values; tol=etol, max_iter=params.maxiters)
        solution, converged = local_newton(params_solver, initial_values; tol=etol, max_iter=params.maxiters, debug)

        # if unconverged, extrapolate instead
        if !converged
            solution = [theta0, deltastar0, n0, ues[2]]
            if airfoil
                solution[1] = theta0 * sqrt(xs[2]/xs[1])
                solution[2] = deltastar0 * sqrt(xs[2]/xs[1])
            else
                rlen = (xs[2]-xs[1]) / (10 * theta0)
                solution[1] = (theta0 + rlen * deltastar0) / (1 + rlen)
            end
            if params.verbose; @warn "BL initialization station not converged. Extrapolating values..."; end
        end
    else

        prob = NonlinearSolve.NonlinearProblem(params_solver, initial_values)
        isoutofdomain = IsOutOfDomain(fill(1e-6, 3)) #todo: generalize better or make a user input?
        sol = NonlinearSolve.solve(prob, NonlinearSolve.LevenbergMarquardt(), 
                                        abstol=etol, 
                                        maxiters=params.maxiters,
                                        isoutofdomain=isoutofdomain)

        if sol.retcode != NonlinearSolve.ReturnCode.Success
            sol2 = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(), 
                                        abstol=etol, 
                                        maxiters=params.maxiters,
                                        isoutofdomain=isoutofdomain)
            if sol2.retcode == NonlinearSolve.ReturnCode.Success
                sol = sol2
            end
        end

        solution = sol.u
    end

    theta_2 = solution[1]
    delta_2 = solution[2]
    n_2 = solution[3]
    ue_2 = solution[4]

    return delta_2, theta_2, n_2, ue_2
end

function station_stagnation(deltastar0, theta0, n0, ues, xs, params; upper=true)

    (; etol, verbose) = params
    
    initial_values = vcat(theta0, deltastar0, n0)

    if any(isnan.(initial_values))
        @infiltrate
    end

    turbulent = false
    transition = false
    similar = true
    airfoil = true #always true for stagnation
    forced_transition_flag = transition
    wake_gap = zeros(2)

    params_solver = ResidualParamsStation(xs, deltastar0, theta0, n0, ues, wake_gap, params, upper, turbulent, transition, similar, airfoil, forced_transition_flag)
    if params.use_native_solver_local
        # solution, converged = local_levenberg_marquardt(params_solver, initial_values; tol=etol, max_iter=params.maxiters)
        #TODO: keep investigating here
        solution, converged = local_newton(params_solver, initial_values; tol=etol, max_iter=params.maxiters)
        # if unconverged, extrapolate instead
        if !converged
            solution = [deltastar0, theta0, n0]
            if airfoil
                solution[1] = theta0 * sqrt(xs[2]/xs[1])
                solution[2] = deltastar0 * sqrt(xs[2]/xs[1])
            else
                rlen = (xs[2]-xs[1]) / (10 * theta0)
                solution[1] = (theta0 + rlen * deltastar0) / (1 + rlen)
            end
            if params.verbose; @warn "BL initialization station not converged. Extrapolating values..."; end
        end
    else

        prob = NonlinearSolve.NonlinearProblem(params_solver, initial_values)
        isoutofdomain = IsOutOfDomain(fill(1e-6, 3)) #todo: generalize better or make a user input?
        sol = NonlinearSolve.solve(prob, NonlinearSolve.LevenbergMarquardt(),
                                    abstol=etol, 
                                    maxiters=params.maxiters,
                                    isoutofdomain=isoutofdomain)

        if sol.retcode != NonlinearSolve.ReturnCode.Success
            sol2 = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(), 
                                        abstol=etol, 
                                        maxiters=params.maxiters,
                                        isoutofdomain=isoutofdomain)
            if sol2.retcode == NonlinearSolve.ReturnCode.Success
                sol = sol2
            end
        end                      

        solution = sol.u
    end

    return solution
end

function station_transition(deltastar1, theta1, n1, ues, xs, params; upper=true, debug=false, forced_transition_flag=false)

    ct2_0 = get_ct_transition(deltastar1, theta1, n1, ues[2], params; turbulent=false)

    # if params.transition_method != 3
        initial_values = vcat(theta1, deltastar1, ct2_0, ues[2])
    # else
    #     initial_values = vcat(theta1, deltastar1, ct2_0, ues[2], (xs[1]+xs[2])/2, theta1, deltastar1)
    # end

    turbulent = false
    transition = true
    similar = false
    airfoil = true
    wake_gap = zeros(2)
    params_solver = ResidualParamsStation(xs, deltastar1, theta1, n1, ues, wake_gap, params, upper, turbulent, transition, similar, airfoil, forced_transition_flag)

    if params.use_native_solver_local
        p = (; params, xs)

        ### ImplictAD implementation - would need a way to tell if it converged or not
        # residuals_transition_implicit(y, x, p) = residuals_transition(p.xs, [x[2], y[2]], [x[1], y[1]], [x[3], y[3]], [x[4], y[4]], p.params; upper, debug=false)
        # local_newton_implicit_temp(x, p) = local_newton_implicit(params_solver, x)
        # solution_implicit = ImplicitAD.implicit(local_newton_implicit_temp, residuals_transition_implicit, initial_values)
        
        # solution, converged = local_levenberg_marquardt(params_solver, initial_values; tol=params.etol, max_iter=params.maxiters)
        solution, converged = local_newton(params_solver, initial_values; debug)

        # if unconverged, extrapolate instead
        if !converged #TODO: make separate function
            solution = initial_values

            solution[1] = theta1 * sqrt(xs[2]/xs[1])
            solution[2] = deltastar1 * sqrt(xs[2]/xs[1])
            if params.verbose; @warn "BL initialization station not converged. Extrapolating values..."; end
        end
    else

        prob = NonlinearSolve.NonlinearProblem(params_solver, initial_values)

        isoutofdomain = IsOutOfDomain(fill(1e-6, 3)) #todo: generalize better or make a user input?
        sol = NonlinearSolve.solve(prob, NonlinearSolve.LevenbergMarquardt(), 
                                        abstol=params.etol, 
                                        maxiters=params.maxiters,
                                        isoutofdomain=isoutofdomain)

        if sol.retcode != NonlinearSolve.ReturnCode.Success
            sol2 = NonlinearSolve.solve(prob, NonlinearSolve.NewtonRaphson(), 
                                        abstol=params.etol, 
                                        maxiters=params.maxiters,
                                        isoutofdomain=isoutofdomain)
            if sol2.retcode == NonlinearSolve.ReturnCode.Success
                sol = sol2
            end
        end

        solution = sol.u
    end

    theta2 = solution[1]
    deltastar2 = solution[2]
    ct2 = solution[3]
    ue2 = solution[4]

    # if params.transition_method == 3
    #     xt = solution[5]
    #     theta_t = solution[6]
    #     delta_t = solution[7]

    #     if upper
    #         params.transition_upper[1] = xt
    #     else
    #         params.transition_lower[1] = xt
    #     end
    # else
        xt = theta_t = delta_t = 0.0
        # todo: at least define these based on the saved xt in params and then use weighting to get the others. 
    # end

    return deltastar2, theta2, ct2, ue2, xt, theta_t, delta_t
end

############ Residual equations ############

function residuals_transition(xs, deltas, thetas, state3s, ues, params; upper=true, debug=false, forced_transition_flag=false)
    (; Mach, etol, ncrit) = params

    U1 = vcat(thetas[1], deltas[1], state3s[1], ues[1])
    U2 = vcat(thetas[2], deltas[2], state3s[2], ues[2])

    x1 = xs[1]
    dx = xs[2] - xs[1]
    
    p = (; params, debug)

    Us = vcat(U1,U2,x1,dx)

    if forced_transition_flag
        xt = upper ? params.sft_u[1] : params.sft_l[1]
    else
        xt = ImplicitAD.implicit(solve_for_xt_implicit, xt_residual_implicit, Us, p)
    end

    if eltype(xt) == Float64
        if upper
            params.transition_upper[1] = xt # set new transition point #todo: for now it will do this multiple times. Will want to improve - maybe based on the norm of Rlam + Rturb and params.etol?
        else
            params.transition_lower[1] = xt
        end
    else
        if upper
            params.transition_upper[1] = xt.value 
        else
            params.transition_lower[1] = xt.value
        end
    end

    #define Ut
    w_turb = (xt - x1) / dx # turbulent weight to apply to states
    w_lam = 1.0 - w_turb # laminar weight to apply to states
    Ut = w_lam * U1 .+ w_turb * U2
    Utl = copy(Ut)
    Utt = copy(Ut)
    if !forced_transition_flag
        Utl[3] = ncrit
    end

    Utt[3] = get_ct_transition(Utl, params; turbulent=true)

    # define new states
    x_lam = [xs[1], xt]
    x_turb = [xt, xs[2]]

    # solve for laminar and turbulent residuals
    Rlaminar = residuals_together_local(x_lam, U1, Utl, params; airfoil=true, turbulent=false, upper)
    Rturbulent = residuals_together_local(x_turb, Utt, U2, params; airfoil=true, turbulent=true, upper)

    if forced_transition_flag
        Rlaminar[3] = 0
    end

    R = Rlaminar .+ Rturbulent

    return R
end

# """
# Solves xt, delta_t, and theta_t in a nested solve before using those solutions to solve for the turbulent states in the global solve
# """
# function residuals_transition_test2(xs, deltas, thetas, state3s, ues, params; upper=true, debug=false)

#  (; Mach, etol, ncrit) = params
#     Me = Mach

#     U1 = vcat(thetas[1], deltas[1], state3s[1], ues[1])
#     U2 = vcat(thetas[2], deltas[2], state3s[2], ues[2])

#     x1 = xs[1]
#     dx = xs[2] - xs[1]
    
#     p = (; params, upper)

#     Us = vcat(U1,U2,x1,dx)

#     # xt = solve_for_xt(Us, p)
#     # xt_implicit = ImplicitAD.implicit(solve_for_xt_implicit, xt_residual_implicit, Us, p)

#     y = ImplicitAD.implicit(transition_nestedsolve_implicit, nested_transition_residual_implicit, Us, p)

#     xt, theta_t, delta_t = y

#     if eltype(xt) == Float64
#         if upper
#             params.transition_upper[1] = xt # set new transition point #todo: for now it will do this multiple times. Will want to improve - maybe based on the norm of Rlam + Rturb and params.etol?
#         else
#             params.transition_lower[1] = xt
#         end
#     else
#         if upper
#             params.transition_upper[1] = xt.value 
#         else
#             params.transition_lower[1] = xt.value
#         end
#     end

#     w_turb = (xt - x1) / dx # turbulent weight to apply to states
#     w_lam = 1.0 - w_turb # laminar weight to apply to states
#     ue_t = w_lam * U1[4] .+ w_turb * U2[4]

#     Utl = [theta_t, delta_t, ncrit, ue_t]

#     ct_t0 = get_ct_transition(Utl, params; turbulent=true)

#     # define new states
#     x_turb = [xt, xs[2]]
#     thetas_turb = [theta_t, U2[1]]
#     deltas_turb = [delta_t, U2[2]]
#     cts_turb = [ct_t0, U2[3]]
#     ues_turb = [ue_t, U2[4]]

#     # solve for turbulent residuals
#     Rturbulent = residuals_together_local(x_turb, deltas_turb, thetas_turb, cts_turb, ues_turb, params; airfoil=true, turbulent=true, upper)
    
#     return Rturbulent
# end

# """
# Solves xt, delta_t, and theta_t as states in the global solve. xs, thetas, and deltas are now length 3 - transition states are the third entries there
# """
# function residuals_transition_test3(xs, deltas, thetas, state3s, ues, params; upper=true, debug=false)

#     (; Mach, etol, ncrit) = params
#     Me = Mach

#     U1 = vcat(thetas[1], deltas[1], state3s[1], ues[1])
#     U2 = vcat(thetas[2], deltas[2], state3s[2], ues[2])

#     x1 = xs[1]
#     dx = xs[2] - xs[1]
#     xt = xs[3]
#     theta_t = thetas[3]
#     delta_t = deltas[3]

#     w_turb = (xt - x1) / dx # turbulent weight to apply to velocity at transition
#     w_lam = 1.0 - w_turb # laminar weight to apply to velocity at transition
#     ue_t = w_lam * U1[4] .+ w_turb * U2[4]

#     Utl = [theta_t, delta_t, ncrit, ue_t]

#     ct_t0 = get_ct_transition(Utl, params; turbulent=true)

#     # define new states
#     x_lam = [xs[1], xt]
#     thetas_lam = [U1[1], theta_t]
#     deltas_lam = [U1[2], delta_t]
#     cts_lam = [U1[3], ncrit]
#     ues_lam = [U1[4], ue_t]
#     x_turb = [xt, xs[2]]
#     thetas_turb = [theta_t, U2[1]]
#     deltas_turb = [delta_t, U2[2]]
#     cts_turb = [ct_t0, U2[3]]
#     ues_turb = [ue_t, U2[4]]

#     # solve for turbulent residuals
#     Rlaminar = residuals_together_local(x_lam, deltas_lam, thetas_lam, cts_lam, ues_lam, params; airfoil=true, turbulent=false, upper)
#     Rturbulent = residuals_together_local(x_turb, deltas_turb, thetas_turb, cts_turb, ues_turb, params; airfoil=true, turbulent=true, upper)
    
#     return vcat(Rlaminar, Rturbulent)
# end

function solve_for_xt(Us, p)

    U1 = Us[1:4]
    U2 = Us[5:8]

    (; x1, dx, params, debug) = p
    p2 = (; x1, dx, params, debug, U1, U2)

    xt, info = FLOWMath.brent(xt_residual, p.x1, p.x1+p.dx; args=(;p2), atol=1e-12)

    if info.flag != "CONVERGED"
        if p.params.verbose; @warn "Transition location solver not converged - splitting it halfway"; end
        xt = x1 + 0.5 * p.dx
    end

    return xt
end

function xt_residual(states, p)

    (; x1, dx, U1, U2, params, debug) = p

    xt = states[1]
    w_turb = (xt - x1) / dx # turbulent weight to apply to states
    w_lam = 1.0 - w_turb # laminar weight to apply to states

    Ut = w_lam .* U1 .+ w_turb .* U2
    
    T = eltype(Ut)
    Ut = setindex!(Ut, T(params.ncrit), 3)  # Explicit type conversion

    dndx1 = get_dndx(U1, params) 
    dndxt = get_dndx(Ut, params) 

    R = params.ncrit - U1[3] - (xt-x1) * (dndx1 + dndxt) / 2

    return R
end

function solve_for_xt_implicit(Us, p)

    U1 = Us[1:4]
    U2 = Us[5:8]
    x1 = Us[9]
    dx = Us[10]

    (; params, debug) = p
    p2 = (; x1, dx, params, debug, U1, U2)

    xt, info = FLOWMath.brent(xt_residual, x1, x1+dx; args=(;p2), atol=1e-12)

    if info.flag != "CONVERGED"
        xt0 = x1 + 0.5*dx
        xt, converged = local_newton_xt(xt_residual, xt0, p2)

        if converged
            if p.params.verbose
                @warn "Transition ($xt) off interval ($x1,$(x1+dx))!"
            end
        else
            if p.params.verbose
                @warn "Transition location solver not converged - splitting it halfway"
            end

            xt = xt0
        end
    end

    return xt
end

function xt_residual_implicit(y, Us, p)

    U1 = Us[1:4]
    U2 = Us[5:8]
    x1 = Us[9]
    dx = Us[10]

    (; params, debug) = p

    xt = y[1]
    w_turb = (xt - x1[1]) / dx # turbulent weight to apply to states
    w_lam = 1.0 - w_turb # laminar weight to apply to states

    Ut = w_lam .* U1 .+ w_turb .* U2
    
    T = eltype(Ut)
    Ut = setindex!(Ut, T(params.ncrit), 3)  # Explicit type conversion

    dndx1 = get_dndx(U1, params) 
    dndxt = get_dndx(Ut, params) 

    return params.ncrit - U1[3] - (xt-x1[1]) * (dndx1 + dndxt) / 2
end

function xt_residual(states, x1, dx, U1, U2, params, debug)

    p = (; x1, dx, U1, U2, params, debug)

    return xt_residual(states, p)
end

residuals_together_local(xs, U1, U2, args...; kwargs...) = residuals_together_local(xs, [U1[2], U2[2]], [U1[1], U2[1]], [U1[3], U2[3]], [U1[4], U2[4]], args...; kwargs...)

function residuals_together_local(xs, deltastars, thetas, state3s, ues, params, wake_gap=zeros(2); airfoil=true, turbulent=true, upper=true, similar=false)
    (; rho, mu, Mach) = params

    if !airfoil #only in the wake
        deltastars .-= wake_gap
    end

    # wake gap parameter - zero for airfoil stations
    Hwake1 = wake_gap[1] / thetas[1]
    Hwake2 = wake_gap[2] / thetas[2]
    Hwake = (Hwake1 + Hwake2) / 2
    
    H_vec = get_H.(deltastars, thetas)

    # determine upwinding average parameter
    M2 = get_Mach2(ues, params)

    Reθ_vec = get_Reθ.(thetas, ues, Ref(params))

    ηup = get_ηup(H_vec[1], H_vec[2], M2[1]; airfoil, turbulent)

    Hk_vec = get_Hk.(H_vec, M2; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)

    Hstarstar_vec = get_Hstarstar.(Hk_vec, M2)

    # mid parameters
    x_mid = xs[1]/2 + xs[2]/2
    δ_mid = deltastars[1]/2 + deltastars[2]/2
    θ_mid = thetas[1]/2 + thetas[2]/2
    H_mid = δ_mid / θ_mid
    ue_mid = ues[1]/2 + ues[2]/2
    M2_mid = get_Mach2(ue_mid, params)

    Hk_mid = get_Hk(H_mid, M2_mid; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)
    Reθ_mid = get_Reθ(θ_mid, ue_mid, params)

    if turbulent
        cτ12s = state3s

        if !airfoil #no cf in the wake
            cf_vec = zeros(2)
            cf_mid = 0.0
        else
            cf_vec = get_cf_turbulent.(Hk_vec, Reθ_vec, M2)
            cf_mid = get_cf_turbulent(Hk_mid, Reθ_mid, M2_mid)
        end
        cfxt_vec = cf_vec .* xs ./ thetas # cf * x / theta
        cfxt_mid = cf_mid * x_mid / θ_mid

        Hstar_vec = get_Hstar_turbulent.(Hk_vec, Reθ_vec, M2; airfoil)

        Us_vec = get_Us.(H_vec, Hk_vec, Hstar_vec; airfoil)

        cτeq12_vec = get_cτ12eq.(H_vec, Hk_vec, Hstar_vec, Reθ_vec, Us_vec; airfoil)

        cdi_vec = get_cdi_turbulent.(Hk_vec, Reθ_vec, cf_vec, Us_vec, Hstar_vec, cτ12s .^ 2; airfoil)
        cdixt_vec = cdi_vec .* xs ./ (thetas)
        # cdixtHstar_vec = cdi_vec .* xs ./ (thetas .* Hstar_vec) #todo: know why we don't include Hstar

        δ_vec = get_δ.(deltastars, thetas, Hk_vec)

        return residuals_together_turbulent(xs, ηup, H_vec, ues, deltastars, thetas, cfxt_vec, cfxt_mid, 
                                            Hstar_vec, Hstarstar_vec,
                                            cdixt_vec, δ_vec, cτ12s, cf_vec, Us_vec, 
                                            cτeq12_vec, Hk_vec, Reθ_vec, Hwake, params; airfoil)
    else
        ñ_vec = state3s

        cf_vec = get_cf_laminar.(Hk_vec, Reθ_vec)
        cf_mid = get_cf_laminar(Hk_mid, Reθ_mid)
        cfxt_vec = cf_vec .* xs ./ thetas # cf * x / theta

        cfxt_mid = cf_mid * x_mid / θ_mid

        Hstar_vec = get_Hstar_laminar.(Hk_vec)

        cdi_vec = get_cdi_laminar.(Hk_vec, Reθ_vec)

        cdixt_vec = cdi_vec .* xs ./ (thetas)# .* Hstar_vec) 
        # cdixt_vec = cdi_vec .* xs ./ (thetas .* Hstar_vec) #todo: come up with a really good explanation for why this should be omitted. Ask Drela and Fidkowski if needed

        R = residuals_together_laminar(xs, ηup, H_vec, ues, thetas, cfxt_vec, cfxt_mid,
                                            Hstar_vec, Hstarstar_vec, 
                                            cdixt_vec, ñ_vec, Hk_vec, Reθ_vec, Hwake, params; airfoil, similar)

        return R
    end
end

function residuals_together_turbulent(xs, ηup, H_vec, ues, deltastars,
                            thetas, cfxt_vec, cfxt_mids,
                            Hstar_vec, Hstarstar_vec, 
                            cdixt_vec, δ_vec, cτ12_vec, cf_vec, Us_vec,
                            cτeq12_vec, Hk_vec, Reθ_vec, Hwake, params; airfoil=true)

    R_mom = residual_momentum(xs, ηup, H_vec, ues, thetas, cfxt_vec, cfxt_mids, Hwake, params; airfoil)

    R_shape = residual_shape(xs, ηup, H_vec, ues, thetas, cfxt_vec, Hstar_vec, 
                            Hstarstar_vec, cdixt_vec, Hwake, params; airfoil)

    R_lag = residual_shearlag(xs, ηup, H_vec, ues, deltastars, δ_vec, cτ12_vec, 
                                Us_vec, cτeq12_vec, cf_vec, Hk_vec, Reθ_vec, params; airfoil)

    return vcat(R_mom, R_shape, R_lag)
end

function residuals_together_laminar(xs, ηup, H_vec, ues, thetas, 
                            cfxt_vec, cfxt_mid, Hstar_vec, Hstarstar_vec, 
                            cdixt_vec, ñ_vec, Hk_vec, Reθ_vec, Hwake, params; airfoil=true, similar=false)
    
    R_mom = residual_momentum(xs, ηup, H_vec, ues, thetas, cfxt_vec, cfxt_mid, Hwake, params; airfoil, similar)

    R_shape = residual_shape(xs, ηup, H_vec, ues, thetas, cfxt_vec, Hstar_vec, 
                            Hstarstar_vec, cdixt_vec, Hwake, params; airfoil, similar)

    R_amp = residual_amplification(xs, H_vec, ues, thetas, ñ_vec, 
                                    Hk_vec, Reθ_vec, params.ncrit; similar)

    return vcat(R_mom, R_shape, R_amp)
end

function residual_momentum(xs, ηup, H_vec, ues, thetas, cfxt_vec, cfxt_mids, Hwake, params; airfoil=true, similar=false)
    M2s = get_Mach2(ues, params)
    if params.Mach > 0
        M2 = apply_upwinding(0.5, M2s[1], M2s[2])
    else
        M2 = M2s
    end

    #apply upwinding
    H = apply_upwinding(0.5, H_vec[1], H_vec[2])

    # uses the central average like XFOIL and MFOIL for better drag accuracy
    cfxt_avg = 0.25 * cfxt_vec[1] + 0.25 * cfxt_vec[2] + 0.5 * cfxt_mids

    uks = get_uk(ues, params)

    if !similar
        logθ = log(thetas[2] / thetas[1])
        logue = log(uks[2] / uks[1])
        logx = log(xs[2] / xs[1])
    else #stagnation case
        logθ = 0.0
        logue = 1.0
        logx = 1.0
    end

    R_mom = logθ + logue * (H + 2 + Hwake - M2) - logx * cfxt_avg / 2
    return R_mom
end

function residual_shape(xs, ηup, H_vec, ues, θ_vec, cfxt_vec, Hstar_vec, Hstarstar_vec, cdixt_vec, Hwake, params; airfoil=true, similar=false)
    uks = get_uk(ues, params)
    #apply upwinding
    H = apply_upwinding(0.5, H_vec[1], H_vec[2]) #todo: determine why Drela uses 0.5, seems arbitrary
    Hstar = apply_upwinding(0.5, Hstar_vec[1], Hstar_vec[2])
    Hstarstar = apply_upwinding(ηup, Hstarstar_vec[1], Hstarstar_vec[2])
    cfxt = apply_upwinding(ηup, cfxt_vec[1], cfxt_vec[2])
    cdixt = apply_upwinding(ηup, cdixt_vec[1], cdixt_vec[2])

    if !similar
        logHstar = log(Hstar_vec[2] / Hstar_vec[1])
        logue = log(uks[2] / uks[1])
        logx = log(xs[2] / xs[1])
    else #stagnation case
        logHstar = 0.0
        logue = 1.0
        logx = 1.0
    end

    R_shape = logHstar + logue * (2*Hstarstar / Hstar + 1 - H - Hwake) + logx * (cfxt/2 - cdixt)
    return R_shape
end

function residual_amplification(xs, H_vec, ues, thetas, ns, Hk_vec, Reθ_vec, ncrit; similar=false)

    if !similar
        dndx_vec = get_dndx.(Hk_vec, Reθ_vec, ns, thetas; ncrit=ncrit)

        dndx = apply_upwinding(0.5, dndx_vec[1], dndx_vec[2])

        R_amp = @. ns[2] - ns[1] - dndx * (xs[2] - xs[1])
    else
        R_amp = ns[2] + ns[1] #stagnation case
    end
    return R_amp
end

function residual_amplification(xs, U1, U2, params)

    dndx1 = get_dndx(U1, params)
    dndx2 = get_dndx(U2, params)

    dndx = apply_upwinding(0.5, dndx1, dndx2)

    R_amp = U2[3] - U1[3] - dndx * (xs[2] - xs[1])
    return R_amp
end

function residual_shearlag(xs, ηup, H_vec, ues, deltastars, δ_vec, cτ12_vec, Us_vec, cτeq12_vec, cf_vec, Hk_vec, Reθ_vec, params; Klag=5.6, GB=0.75, airfoil=true)
    uks = get_uk(ues, params)
    ηD = airfoil ? 1.0 : 0.9 # 0.9 in wake

    Δx = xs[2] .- xs[1]

    #apply upwinding
    δstar_avg = deltastars[1]/2 .+ deltastars[2]/2
    δ = apply_upwinding(0.5, δ_vec[1], δ_vec[2])
    Us = apply_upwinding(0.5, Us_vec[1], Us_vec[2])
    cτ12 = apply_upwinding(ηup, cτ12_vec[1], cτ12_vec[2])
    cτeq12 = apply_upwinding(ηup, cτeq12_vec[1], cτeq12_vec[2])
    cf = apply_upwinding(ηup, cf_vec[1], cf_vec[2])
    Hk = apply_upwinding(ηup, Hk_vec[1], Hk_vec[2])
    Reθ = apply_upwinding(0.5, Reθ_vec[1], Reθ_vec[2])

    uq = get_uq(Hk, cf, Reθ, δstar_avg; airfoil)

    R_shearlag = 2δ * log(cτ12_vec[2] / cτ12_vec[1]) - Klag / (GB * (1 + Us)) * Δx * (cτeq12 - ηD * cτ12) - 2δ * (uq * Δx - log(uks[2]/uks[1]))
    return -R_shearlag
end

############ Weighted average with upwinding ############

"""
Apply upwinding weighted average to any given two parameters.
"""
function apply_upwinding(ηup, x1, x2) 
    return (1 - ηup) * x1 + ηup * x2
end

"""
Obtain the upwinding factor for applying a weighted average between panel nodes for certain parameters.
"""
function get_ηup(H1, H2, M2; airfoil=true, turbulent=true)

    Hk1 = get_Hk(H1, M2; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)
    Hk2 = get_Hk(H2, M2; airfoil, turbulent, limit_Hk_max=false, limit_Hk_min=false)

    fHu = (Hk2 - 1) / (Hk1 - 1)

    Cup = airfoil ? 5 : 1 # from xfoil & mfoil code (mfoil paper has a typo)

    return 1 - 0.5 * exp(-(log(abs(fHu)))^2 * Cup/Hk2^2)
end

############ General closure equations ############

function initialize_thermo(Re, Mach, rho, chord, Vinf)
    mu = rho*Vinf*chord/Re
    γair = 1.4      #ratio of specific heats
    γm1 = γair - 1
    KTb = 1.0         #Karman-Tsien beta
    KTl = 0.0         #Karman-Tsien lambda
    H0 = 0.0          #stagnation enthalpy
    mu0 = mu         #stagnation viscosity
    rho0 = rho        #stagnation density
    Tsrat = 0.35    #Suttherland Ts/Tref
    if Mach > 0
        KTb = sqrt(1 - Mach^2)
        KTl = Mach^2/(1 + KTb)^2
        H0 = (1 + 0.5 * γm1 * Mach^2) * Vinf^2 / (γm1 * Mach^2)
        Tr = 1 - 0.5 * Vinf^2 / H0
        finf = Tr^1.5 * (1 + Tsrat)/(Tr + Tsrat)
        mu0 = mu/finf
        rho0 = rho * (1 + 0.5*γm1*Mach^2)^(1/γm1)
    end
    return mu, KTb, KTl, H0, mu0, rho0
end
# get_Me(ue; a=100.0) = (ue / a)
get_Me(ue; a=343.0) = (ue / a)
"""
Useful for compressible, non-zero Mach
"""
function get_uk(ue, params)
    (; Vinf, Mach, KTl) = params

    if Mach > 0
        den = 1 .- KTl .* (ue/Vinf).^2
        uk = ue * (1-KTl) ./ den
    else
        uk = ue
    end

    return uk
end

function get_Mach2(ue, params)  #::OperatingParameters
    (;Mach, γair, H0) = params
    uk = get_uk(ue, params)
    if Mach > 0
        γmi = γair - 1
        c2 = γmi .* (H0 .- 0.5 .* uk .^ 2)
        Mach2 = uk .^ 2 ./c2
    else
        Mach2 = 0
    end

    return Mach2
end

function get_rho(ues, params)
    Mach = params.Mach
    γmi = params.γair - 1
    uks = get_uk(ues, params)
    if Mach > 0
        M2 = get_Mach2(uks, params)
        den = 1 + 0.5 * γmi * M2
        rho = params.rho0/den^(1/γmi)
    else
        rho = param.rho0
    end
    return rho
end

get_Reθ(ρ, ue, θ, μ) = ρ * ue * θ / μ
function get_Reθ(θ, ue, params)
    (; Mach, rho, mu, H0, γair, rho0, mu0) = params
    Ts = 0.35           #sutherland's ratio
    uk = get_uk(ue, params)
    if Mach > 0
        M2 = get_Mach2(ue, params)
        Tr = 1 - 0.5 * uk^2/H0
        f = Tr^1.5 * (1 + Ts)/(Tr + Ts)
        mu1 = mu0*f
        den = 1 + 0.5 * (γair - 1) * M2
        rho1 = rho0 / den^(1/(γair - 1))
        Reθ = get_Reθ(rho1, uk, θ, mu1)
    else
        Reθ = get_Reθ(rho, uk, θ, mu)
    end
    return Reθ
end

function get_H(delta::T, theta::T) where {T<:Number}
    return delta / theta
end

# Both are Duals
function get_H(delta::ForwardDiff.Dual{Tag,T,N},
               theta::ForwardDiff.Dual{Tag,T,N}) where {Tag,T,N}

    H_val = delta.value / theta.value
    dH_ddelta = 1 / theta.value
    dH_dtheta = -H_val / theta.value

    partials = ntuple(i -> dH_ddelta * delta.partials[i] +
                         dH_dtheta * theta.partials[i], N)

    return ForwardDiff.Dual{Tag}(H_val, partials)
end

# delta is Dual, theta is Real
function get_H(delta::ForwardDiff.Dual{Tag,T,N},
               theta::T) where {Tag,T,N}
    H_val = delta.value / theta
    dH_ddelta = 1 / theta

    partials = ntuple(i -> dH_ddelta * delta.partials[i], N)

    return ForwardDiff.Dual{Tag}(H_val, partials)
end

# delta is Real, theta is Dual
function get_H(delta::T,
               theta::ForwardDiff.Dual{Tag,T,N}) where {Tag,T,N}
    H_val = delta / theta.value
    dH_dtheta = -H_val / theta.value

    partials = ntuple(i -> dH_dtheta * theta.partials[i], N)

    return ForwardDiff.Dual{Tag}(H_val, partials)
end

function get_Hk(H, M2; airfoil=true, turbulent=true, limit_Hk_max=true, limit_Hk_min=true)

    Hk = (H .- 0.29 .* M2) / (1 .+ 0.113 .* M2)

    Hk = restrict_Hk(Hk; airfoil, turbulent, limit_Hk_max, limit_Hk_min)

    return Hk
end

function get_Hk(δ, θ, params::OperatingParameters; kwargs...)
    M2 = (params.Mach)^2
    H = δ / θ
    return get_Hk(H, M2; kwargs...)
end

function get_Hk(δ, θ, ue, params::OperatingParameters; kwargs...)
    M2 = get_Mach2(ue, params)
    H = δ / θ
    return get_Hk(H, M2; kwargs...)
end

function restrict_Hk(Hk; airfoil=true, turbulent=true, limit_Hk_min=false, limit_Hk_max=false)

    if limit_Hk_min
        if airfoil
            Hk = max(Hk, 1.05)
        else
            Hk = max(Hk, 1.00005)
        end
    end

    if limit_Hk_max
        if turbulent
            if Hk > 2.5
                Hk = 2.5
                separation_flag = true
            end
        else
            if Hk > 3.8
                Hk = 3.8
                separation_flag = true
            end
        end
    end

    return Hk
end

"""
BL thickness, as defined by mfoil code (not defined in the paper)
"""
function get_δ(deltastar, theta, Hk)
    return min(12*theta, deltastar + theta * (3.15 + 1.72 / (Hk-1)))
end

"""
H**: Density shape parameter
"""
function get_Hstarstar(Hk, M2)
    return M2 * (0.064/ (Hk - 0.8) + 0.251)
end

"""
H*_laminar: Kinetic energy shape parameter

This is by definition δKE / θ
"""
function get_Hstar_laminar(Hk)

    Hk = restrict_Hk(Hk; airfoil=true, limit_Hk_min=true)

    H̃k = Hk - 4.35

    if Hk < 4.35
        Hstar = (0.0111 * H̃k^2 - 0.0278 * H̃k^3) / (Hk + 1) + 1.528 - 0.0002 * (H̃k * Hk)^2
    else
        Hstar = 0.015 * (H̃k)^2 / Hk + 1.528
    end

    return Hstar
end

function get_Hstar_laminar(δstar, θ, ñ, ue, params)

    H = δstar/θ

    M2 = get_Mach2(ue, params)

    Hk = get_Hk(H, M2)

    return get_Hstar_laminar(Hk)
end

get_Hstar_laminar(U, params) = get_Hstar_laminar(U[2], U[1], U[3], U[4], params)


"""
H*_turbulent: Kinetic energy shape parameter

This is by definition δKE / θ
"""
function get_Hstar_turbulent(Hk, Reθ, M2; airfoil=true)

    Hk = restrict_Hk(Hk; airfoil, limit_Hk_min=true)

    Reθ_tilde = max(Reθ, 200) #todo: smooth
    H0 = min(3 + 400/Reθ, 4)
    AHstar = Hk - H0 + 4/log(Reθ_tilde)
    Hr = (H0 - Hk) / (H0 - 1)

    if Hk < H0
        Hstar_turb_incomp = 1.5 + 4/Reθ_tilde + 1.5 * Hr^2 / (Hk + 0.5) * (0.5 - 4/Reθ_tilde)
    else
        Hstar_turb_incomp = 1.5 + 4/Reθ_tilde + (Hk - H0)^2 * (0.015/Hk + 0.007*log(Reθ_tilde)/AHstar^2)
    end

    Hstar = (Hstar_turb_incomp + 0.028 * M2) / (1 + 0.014 * M2)

    return Hstar
end

function get_Hstar_turbulent(δstar, θ, cτ, ue, params; kwargs...)

    M2 = get_Mach2(ue, params)
    H = δstar/θ
    Hk = get_Hk(H, M2)
    Reθ = get_Reθ(θ, ue, params)

    return get_Hstar_turbulent(Hk, Reθ, M2; kwargs...)
end

# get_Hstar_turbulent(U, Me=0.0; kwargs...) = get_Hstar_turbulent(U[2], U[1], U[3], U[4], Me; kwargs...)

"""
cf_laminar: friction coefficient in laminar flow
"""
function get_cf_laminar(Hk, Reθ)

    if real(Hk) < 5.5
        cf_lam = 0.0727 * (5.5 - Hk)^3 / (Hk + 1) - 0.07
    else
        cf_lam = 0.015 * (1 - 1/(Hk - 4.5))^2 - 0.07
    end

    return cf_lam / Reθ
end

"""
cf_wake: friction coefficient in the wake
"""
get_cf_wake() = 0

"""
cf_turbulent: friction coefficient in turbulent flow
"""
function get_cf_turbulent(Hk, Reθ, M2; γair=1.4)

    Fc = sqrt(1 + 0.5 * (γair - 1) * M2)

    Acf = -1.33 * Hk

    if real(Acf) < -17
        Acf = -20 + 3 * exp((Acf+17)/3)
    end

    Bcf = max(log10(Reθ/Fc), 1.303)

    cf = 0.3 * exp(Acf) * Bcf^(-1.74-0.31*Hk) + 0.00011 * (tanh(4 - Hk/0.875) - 1)

    return cf / Fc
end

"""
cdi_laminar: dissipation coefficient in laminar flow
"""
function get_cdi_laminar(Hk, Reθ)

    if real(Hk) < 4
        cdi_lam = 0.00205 * (4 - Hk)^5.5 + 0.207
    else
        cdi_lam = -(0.0016 * (Hk-4)^2) / (1 + 0.02 * (Hk - 4)^2) + 0.207
    end

    return cdi_lam / Reθ
end

# get_cdi_turbulent(cdi_wall, cdi_outer, cdi_stress, cdi_laminar) = max(cdi_wall + cdi_outer + cdi_stress, cdi_laminar)

function get_cdi_turbulent(Hk, Reθ, cf, Us, Hstar, cτ; airfoil=true) 
    
    cdi_wall = airfoil ? get_cdi_wall(Hk, Reθ, cf, Us, Hstar) : 0 #0 in the wake
    cdi_outer = get_cdi_outer(cτ, Hstar, Us)
    cdi_stress = get_cdi_stress(Reθ, Hstar, Us)
    cdi_laminar = airfoil ? get_cdi_laminar(Hk, Reθ) : get_cdi_lamwake(Hk, Reθ)

    scale = airfoil ? 1.0 : 2.0

    return max(cdi_wall + cdi_outer + cdi_stress, cdi_laminar) * scale
end

"""
Turbulent wall contribution to dissipation coefficient
"""
function get_cdi_wall(Hk, Reθ, cf, Us, Hstar)
    return (cf * Us / (2Hstar)) * (1 + tanh(((Hk - 1) * log(Reθ)) / 2.1))
end

"""
Turbulent outer layer contribution to dissipation coefficient
"""
function get_cdi_outer(cτ, Hstar, Us)
    return 2 * cτ * (0.995 - Us) / Hstar
end

"""
Laminar stress contribution to dissipation coefficient
"""
function get_cdi_stress(Reθ, Hstar, Us)
    return 0.3 * (0.995 - Us)^2 / (Hstar * Reθ)
end
"""
Laminar wake dissipation coefficient
"""
function get_cdi_lamwake(Hk, Reθ)
    Hstar = get_Hstar_laminar(Hk) # force laminar

    return 2.2 * (1 - 1/Hk)^2 / (Hk * Hstar * Reθ)
end

cdi_wake(cdi_outer, cdi_stress, cdi_lamwake) = 2 * max(cdi_outer + cdi_stress, cdi_lamwake)

############ Shear stress lag closure equations ############

"""
Normalized wall slip velocity
"""
function get_Us(H, Hk, Hstar; GB=0.75, airfoil=true)

    Hk = restrict_Hk(Hk; airfoil, limit_Hk_min=true)

    Us = Hstar / 2 * (1 - (Hk - 1) / (H * GB))

    if airfoil
        return min(Us, 0.98)
    else
        return min(Us, 0.99995)
    end
end

function get_uq(Hk, cf, Reθ, δstar; GA=6.7, GB=0.75, airfoil=true)

    Hkc = get_Hkc(Hk, Reθ; airfoil)

    ηD = airfoil ? 1.0 : 0.9

    return (0.5cf - (Hkc / (GA * ηD * Hk))^2) / (GB * δstar)
end

function get_Hkc(Hk, Reθ; GC=18.0, airfoil=true)

    GC = airfoil ? GC : 0.0 # 0 in the wake #TODO: smooth?

    return max(Hk - 1 - GC / Reθ, 0.01)
end

"""
Equilibrium value for cτ (square root value)
"""
function get_cτ12eq(H, Hk, Hstar, Reθ, Us; GA=6.7, GB=0.75, airfoil=true)

    # limit Hk
    Hk = restrict_Hk(Hk; airfoil, limit_Hk_min=true)
    
    Hkc = get_Hkc(Hk, Reθ; airfoil)

    # local cτ12eq
    # try
    # @infiltrate
        cτ12eq = sqrt(Hstar * (Hk - 1) * Hkc^2 / (2 * GA^2 * GB * (1 - Us) * H * Hk^2))
    # catch
    #     cτ12eq = sqrt(Complex(Hstar * (Hk - 1) * Hkc^2 / (2 * GA^2 * GB * (1 - Us) * H * Hk^2)))
    # end

    return cτ12eq
end

function get_ct_transition(delta, theta, n, ue, params; turbulent=false)

    (; Cτ0, Eτ0) = params

    M2 = get_Mach2(ue, params)
    H = delta / theta
    Hk = get_Hk(H, M2; turbulent, limit_Hk_max=false)
    Reθ = get_Reθ(theta, ue, params)

    if turbulent
        Hstar = get_Hstar_turbulent(Hk, Reθ, params.Mach; airfoil=true)
    else
        Hstar = get_Hstar_laminar(Hk)
    end

    Us = get_Us(H, Hk, Hstar)

    local cτeq
    try
        cτeq = get_cτ12eq(H, Hk, Hstar, Reθ, Us)
    catch
        @infiltrate
    end
    ct_transition = Cτ0 * exp(-Eτ0 / (Hk - 1)) * cτeq

    return ct_transition
end

get_ct_transition(U, params; turbulent=false) = get_ct_transition(U[2], U[1], U[3], U[4], params; turbulent)

############ Laminar amplification closure equations ############
function get_dndx(rn, fn, gn, ϵn, θ)
    return (rn * fn * gn + ϵn) / θ
end

function get_dndx(Hk, Reθ, n, θ; ncrit=9)

    Hk = restrict_Hk(Hk; airfoil=true, limit_Hk_min=true)

    Ĥ = get_Ĥ(Hk)
    L0 = get_L0(Ĥ)
    sn = get_sn(Reθ, L0)
    rn = get_rn(sn)
    fn = get_fn(Ĥ)
    gn = get_gn(Ĥ)
    ϵn = get_ϵn(n, ncrit)

    return get_dndx(rn, fn, gn, ϵn, θ)
end

function get_dndx(U, params)

    M2 = get_Mach2(U[4], params)
    H = U[2] / U[1]
    Hk = get_Hk(H, M2; turbulent=false, limit_Hk_min=true, limit_Hk_max=false)
    Reθ = get_Reθ(U[1], U[4], params)

    return get_dndx(Hk, Reθ, U[3], U[1]; ncrit=params.ncrit)
end

function get_fn(Ĥ)
    return -0.05 + 2.7Ĥ - 5.5Ĥ^2 + 3Ĥ^3 + 0.1 * exp(-20Ĥ)
end

function get_gn(Ĥ)
    return 0.028 / Ĥ - 0.0345 * exp(-(3.87 * Ĥ - 2.52)^2)
end

function get_rn(sn)
    if real(sn) < 0
        return 0
    elseif real(sn) >= 1
        return 1
    else
        return 3sn^2 - 2sn^3
    end
end

function get_sn(Reθ, L0)
    sn = (log10(Reθ) - (L0 - 0.1)) / 0.2
    return sn
end

function get_L0(Ĥ)
    return 2.492 * Ĥ^0.43 + 0.7 * (1 + tanh(14 * Ĥ - 9.24))
end

function get_Ĥ(Hk)
    return 1 / (Hk - 1)
end

function get_ϵn(n, ncrit)
    return 0.001 * (1 + tanh(5 * (n - ncrit)))
end

"""
    interpolate_stagnation_state(U1, U2, x1, x2)

Initialize the stagnation point states, which are then used to solve for 
    subsequent states along the boundary layer. Needs to be called each 
    iteration in the global solve.
"""
function interpolate_stagnation_state(U1, U2, x1, x2)

    dx = x2 - x1
    rx = x2/x1
    w1 = x2/dx
    w2 = -x1/dx

    Ust = U1*w1 + U2*w2

    # quadratic extrapolation of edge velocity
    wk1 = rx/dx
    wk2 = -1/(rx*dx)
    K = wk1 * U1[4] + wk2 * U2[4]
    xst = eltype(U1)(1e-6) #ensure proper type
    
    Ust[4] = K*xst

    #TODO: redo this to make it less weird
    if eltype(Ust) == Float64 
        if any(Ust .< 0) # this will trigger sometimes if value is nonnegative but derivatives are negative
            idx = findall(Ust .< 0)

            if 3 in idx
                Ust[3] = 0.0
            else
                # # @infiltrate #haven't seen this case yet #todo: make more robust for other negative indices (for loop)
                # for idx_local in idx
                #     if idx_local == 3
                #         Ust[idx_local] = 0.0
                #     else
                #         Ust[idx_local] *= -1
                #     end
                # end
            end
        end
    else
        if any([z.value for z in Ust] .< 0)
            idx = findall([z.value for z in Ust] .< 0)

            if 3 in idx
                Ust[3] -= Ust[3].value
            else
                # @infiltrate #todo: make more robust for other negative indices (for loop)
            end
        end
    end 

    return Ust, xst
end

"""
    initialize_wake(UTE_l, UTE_u, ue1, panel_TE, is_turbulent_l, is_turbulent_u, operatingparameters)

Return the state at the first wake node based on the states at the trailing edge.
"""
function initialize_wake_first_node(UTE_l, UTE_u, Uw, panel_TE, 
                        TE_l_isturbulent, TE_u_isturbulent, 
                        operatingparameters; init_flag=false)

    turbulent = true

    if TE_l_isturbulent
        ctau_l = UTE_l[3]
    else
        ctau_l = get_ct_transition(UTE_l, operatingparameters; turbulent)
    end

    if TE_u_isturbulent
        ctau_u = UTE_u[3]
    else
        ctau_u = get_ct_transition(UTE_u, operatingparameters; turbulent)
    end

    theta_l = UTE_l[1]
    theta_u = UTE_u[1]
    deltastar_l = UTE_l[2]
    deltastar_u = UTE_u[2]

    deltastar_wake = deltastar_l + deltastar_u + panel_TE.hTE
    theta_wake = theta_l + theta_u
    ctau_wake = (theta_l * ctau_l + theta_u * ctau_u) / (theta_l + theta_u)

    if init_flag # return U
        return [theta_wake, deltastar_wake, ctau_wake]
    else
        R_shape = Uw[2] - deltastar_wake #return R
        R_mom = Uw[1] - theta_wake
        R_lag = Uw[3] - ctau_wake

        return [R_mom, R_shape, R_lag]
    end
end

################ calculating final forces #################

function get_cl(panels_airfoil, opts, params)
    xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[end].x2)
    ys = vcat([p.y1 for p in panels_airfoil], panels_airfoil[end].y2)

    ues = vcat(opts.Vel, opts.Veu)

    return get_cl(xs, ys, ues, params.Vinf, sin(params.alpha[1]), cos(params.alpha[1]))
end

function get_cl(xs, ys, ues, Vinf, sinalpha, cosalpha)
    
    cps = 1 .- (ues ./ Ref(Vinf)) .^ 2
    npanels = length(xs)-1

    cl = 0
    for i in 1:npanels
        cp_bar = (cps[i+1]+cps[i])/2
        dx = xs[i+1]-xs[i]
        dy = ys[i+1]-ys[i]
        ds = -dx*cosalpha - dy*sinalpha

        cl += cp_bar*ds
    end

    cp_bar_end = (cps[1]+cps[npanels+1])/2
    dxend = xs[1]-xs[end]
    dyend = ys[1]-ys[end]
    dsend = -dxend*cosalpha - dyend*sinalpha

    cl += cp_bar_end*dsend

    return cl
end

function get_cps(opts, params)
    (; Vel, Veu, Vwake) = opts
    (; Vinf, Mach, KTb, KTl) = params

    cp_l = 1 .- (Vel ./ Vinf) .^ 2
    cp_u = 1 .- (Veu ./ Vinf) .^ 2
    cp_w = 1 .- (Vwake ./ Vinf) .^ 2
    if (Mach > 0)
        den_l = KTb .+ (0.5*KTl*(1+KTb)) .* cp_l
        cp_l = cp_l ./ den_l

        den_u = KTb .+ (0.5*KTl*(1+KTb)) .* cp_u
        cp_u = cp_u ./ den_u

        den_w = KTb .+ (0.5*KTl*(1+KTb)) .* cp_w
        cp_w = cp_w ./ den_w
    end

    cps = vcat(cp_l, cp_u, cp_w)

    return (; cps, cp_l, cp_u, cp_w)
end

function get_cd(opts, params)
    θ = opts.thetas_w[end]
    ue = opts.Vwake[end]
    H = (opts.deltas_w[end])/θ
    V∞ = params.Vinf
    return get_cd(θ, ue, H, V∞)
end

function get_cd(θ, ue, H, V∞)
    return 2*θ*(ue/V∞)^((5+H)/2)
end

function get_cf_af(opts, params)

    idx_transition_u = FLOWMath.findindex(opts.su, opts.stu)
    idx_transition_l = FLOWMath.findindex(opts.sl, opts.stl)

    cf_u = get_cf_surface(params, opts.deltas_u, opts.thetas_u, opts.Veu, idx_transition_u)
    cf_l = get_cf_surface(params, opts.deltas_l, opts.thetas_l, opts.Vel, idx_transition_l)

    return (; cf_u, cf_l)
end

function get_cf_surface(params, deltas, thetas, ues, idx_transition)
    H = deltas ./ thetas
    Reθ = get_Reθ.(thetas, ues, Ref(params))
    M2 = get_Mach2.(ues, Ref(params))
    
    Hk_laminar = get_Hk.(H[1:idx_transition-1], M2[1:idx_transition-1]; turbulent=false, limit_Hk_max=false)
    Hk_turb = get_Hk.(H[idx_transition:end], M2[idx_transition:end]; turbulent=true, limit_Hk_max=false)
    
    cf_laminar = get_cf_laminar.(Hk_laminar, Reθ[1:idx_transition-1])
    cf_turbulent = get_cf_turbulent.(Hk_turb, Reθ[idx_transition:end], M2[idx_transition:end])

    cf = vcat(cf_laminar, cf_turbulent)
    return cf
end

# function get_cdf_mfoil(deltas, thetas, ues, turbulent, upper, lower, x, y, stag; rho=1, Vinf=1, Mach=0.0, mu=1e-5, chord=1.0, alpha=0.0)

#     idx_transition_u = something(findfirst(x -> x>0, turbulent[upper[1]:upper[2]]), length(turbulent[upper[1]:upper[2]]))
#     idx_transition_l = something(findfirst(x -> x>0, turbulent[lower[1]:-1:lower[2]]), length(turbulent[lower[1]:-1:lower[2]]))

#     cfl = get_cf_surface(rho, mu, Mach, deltas[lower[1]:-1:lower[2]], thetas[lower[1]:-1:lower[2]], ues[lower[1]:-1:lower[2]], idx_transition_l)
#     cfu = get_cf_surface(rho, mu, Mach, deltas[upper[1]:upper[2]], thetas[upper[1]:upper[2]], ues[upper[1]:upper[2]], idx_transition_u)

#     stag_point = stag[1] .* cos(alpha) + stag[2] .* sin(alpha)
#     x_dir = x .* cos(alpha) + y .* sin(alpha)
#     xu = vcat(stag_point, x_dir[upper[1]:upper[2]])
#     xl = vcat(stag_point,x_dir[lower[1]:-1:lower[2]])

#     return get_cdf(cfu, cfl, ues[upper[1]:upper[2]], ues[lower[1]:-1:lower[2]], xu, xl, rho, Vinf, chord)
# end

function get_cdf(cf_u, cf_l, ue_u, ue_l, xu, xl, ρ, V∞, c)
    τw_u = (1/2) * ρ .* cf_u .* ue_u.^2
    D_u = trapz(xu, vcat(0, τw_u))
    τw_l = (1/2) * ρ .* cf_l .* ue_l.^2
    D_l = trapz(xl, vcat(0, τw_l))
    return (D_u + D_l)/((1/2)*ρ*c*V∞^2)
end

function get_cdf(panels_airfoil, opts, params)
    (;rho, Vinf, chord, alpha) = params
    (; s, Veu, Vel, s_stag) = opts

    cf_u, cf_l = get_cf_af(opts, params)

    xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[end].x2)
    ys = vcat([p.y1 for p in panels_airfoil], panels_airfoil[end].y2)

    idx_stag = FLOWMath.linear(s, 1:length(s), s_stag)
    idx_m1 = Int(floor(idx_stag))
    idx_p1 = Int(ceil(idx_stag))
    xs_stag = FLOWMath.linear([idx_m1, idx_p1], [xs[idx_m1], xs[idx_p1]], idx_stag)
    ys_stag = FLOWMath.linear([idx_m1, idx_p1], [ys[idx_m1], ys[idx_p1]], idx_stag)

    xu = vcat(xs_stag, xs[idx_p1:end])
    xl = vcat(xs_stag, xs[idx_m1:-1:1])

    yu = vcat(ys_stag, ys[idx_p1:end])
    yl = vcat(ys_stag, ys[idx_m1:-1:1])

    xu_dir = xu .* cos(alpha[1]) + yu .* sin(alpha[1])
    xl_dir = xl .* cos(alpha[1]) + yl .* sin(alpha[1])

    return get_cdf(cf_u, cf_l, Veu, Vel, xu_dir, xl_dir, rho, Vinf, chord)
end

function get_cdp(panels_airfoil, opts, params)
    cd = get_cd(opts, params)
    cdf = get_cdf(panels_airfoil, opts, params)
    return cd - cdf, cd
end

function initialize_transition_velocity(xt, s, ues)
    # interpolate velocity at transition
    return FLOWMath.linear(s, ues, xt)
end