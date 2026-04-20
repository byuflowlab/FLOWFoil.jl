#todo: pass in wrapper instead of redefining it each time based on input p
function local_newton(F, x0; tol=F.params.etol, max_iter=F.params.maxiters, under_relaxation=true, amp_only=false, debug=false, iters_direct=12)
    x = copy(x0)
    direct = true
    Hk_step = 0.0
    debug = false
    ds = F.xs[2] - F.xs[1]
if any(isnan.(x0))
    @infiltrate
end
    get_Hk_wrapper(states) = get_Hk(states[2], states[1], F.ues[2], F.params; airfoil=F.airfoil, turbulent=F.turbulent, limit_Hk_max=false, limit_Hk_min=false)

    for i in 1:max_iter
        local Fx
        Fx = F(x)

        if norm(Fx) < tol # Convergence check
            return x, true
        end
        if direct
            if length(x) == 1
                J = ForwardDiff.gradient(F, x)
            elseif F.transition && F.params.transition_method==3
                xtemp = vcat(x[1:3], x[5:7])
                J = ForwardDiff.jacobian(F, xtemp)
            else
                xtemp = x[1:3]
                J = ForwardDiff.jacobian(F, xtemp)  # Automatic Jacobian
            end
            local Δx
            try
                Δx = pinv(J) * -Fx  # Linear system solve
            catch
                @infiltrate
            end
        else

            # x4 = vcat(x, F.ues[2])
            J1 = ForwardDiff.jacobian(F, x) 
            J2 = ForwardDiff.gradient(get_Hk_wrapper, x)
            if debug; @infiltrate; end
            if F.transition && F.params.transition_method==3
                Jnew = vcat(J1[1:3,:], J2', J1[4:6,:])
            else
                Jnew = vcat(J1, J2') #should be 4x4
            end

            Hk = get_Hk_wrapper(x)

            if F.transition && F.params.transition_method==3
                # @infiltrate
                # debug = true
                b = vcat(-Fx[1:3], Hk_step - Hk, -Fx[4:6])
            else
                b = vcat(-Fx, Hk_step - Hk)
            end

            Δx = pinv(Jnew) * b
        end

        #under relaxation
        if under_relaxation
            if F.transition && F.params.transition_method==3
                if direct
                    dm = max(abs(Δx[6]/F.δstar0), abs(Δx[5]/F.θ0)) # laminar to transition nodes
                    dm = max(dm, abs(Δx[6]/(x[2]+Δx[2])), abs(Δx[5]/(x[1]+Δx[1]))) # transition to turbulent nodes
                else
                    dm = max(abs(Δx[7]/F.δstar0), abs(Δx[6]/F.θ0)) # laminar to transition nodes
                    dm = max(dm, abs(Δx[7]/(x[2]+Δx[2])), abs(Δx[6]/(x[1]+Δx[1]))) # transition to turbulent nodes
                end
            else # laminar to turbulent nodes
                dm = max(abs(Δx[2]/F.δstar0), abs(Δx[1]/F.θ0))
            end
            if !direct
                dm = max(dm, abs(Δx[4]/F.ues[1]))
            end
            if F.turbulent
                dm = max(dm, abs(Δx[3]/F.state30))
            # elseif F.transition && F.params.transition_method==3
            #     w_turb = ((x[5]+Δx[4]) - F.xs[1]) / ds # turbulent weight to apply to states
            #     w_lam = 1.0 - w_turb # laminar weight to apply to states
            #     ue_t = w_lam * F.ues[1] .+ w_turb * x[4]
            #     Utl_temp = [x[6]+Δx[5], x[7]+Δx[6], F.params.ncrit, ue_t]
            #     ct0 = get_ct_transition(Utl_temp, F.params; turbulent=true)
            #     dm = max(dm, abs(Δx[3]/ct0))
            elseif direct
                dm = max(dm, abs(Δx[3]/10))
            end
            omega = dm > 0.3 ? 0.3/dm : 1.0
if debug; @infiltrate; end
            Δx .*= omega
        end

        if length(x) == 1 
            xnew = x + Δx
        elseif (F.similar || !direct)
            xnew = x + Δx
        elseif F.transition && F.params.transition_method==3
            xnew = vcat(x[1:3] + Δx[1:3], x[4], x[5:7] + Δx[4:6])
        else
            xnew = vcat(x[1:3] + Δx, x[4]) #todo: add new case
        end
if debug; @infiltrate; end
        # Prevent extreme values
        if !amp_only
            if F.turbulent || F.transition
                if xnew[3] > 0.3
                    xnew[3] = 0.3
                elseif xnew[3] < 1e-7
                    xnew[3] = 1e-7
                end
                if F.transition && F.params.transition_method==3
                    # make sure it doesn't jump more than a single panel away?
                    if xnew[5] < (F.xs[1] - ds)
                        xnew[5] = F.xs[1] - ds
                    elseif xnew[5] > (F.xs[2] + ds)
                        xnew[5] = F.xs[2] + ds
                    end

                    # if xnew[5] < F.xs[1] #xt needs to be within the panel
                    #     xnew[5] = F.xs[1]
                    # elseif xnew[5] > F.xs[2]
                    #     xnew[5] = F.xs[2]
                    # end
                end
            else
                if xnew[3] < 0
                    xnew[3] = 0
                end
            end
        else
            if !F.turbulent
                if xnew[1] < 0
                    xnew[1] = 0
                end 
            end
        end
if debug; @infiltrate; end
        ####### check for separation #######
        Hk = get_Hk_wrapper(xnew)
        Hkmax = F.turbulent ? 2.5 : 3.8

        if direct && ((Hk > Hkmax) || i > iters_direct || any(xnew[1:2] .< 0)) 
            #don't take the update; prescribe Hk instead; skip if already separated

            direct = false
            Hk = get_Hk(F.δstar0, F.θ0, F.ues[1], F.params; airfoil=F.airfoil, turbulent=F.turbulent, limit_Hk_max=false, limit_Hk_min=false)
            Hkr = ds / F.θ0

            if !F.airfoil
                H2 = Hk
                for k in 1:6
                    H2 = H2 - (H2 + 0.03 * Hkr * (H2-1)^3 - Hk) / (1 + 0.09 * Hkr * (H2-1)^2)
                end
                Hk_step = max(H2, 1.01)
            else
                Hk_step = F.turbulent ? Hk - 0.15*Hkr : Hk + 0.03 * Hkr
            end
            if F.airfoil
                Hk_step = max(Hk_step, Hkmax)
            end

            # if i > iters_direct && !debug
            #     x = copy(x0) #reinitialize states
            # end

        else
            x = xnew #take the update
        end
    end

    # check one more time
    Fx = F(x)

    if norm(Fx) < tol # Convergence check
        return x, true
    else
        return x, false
    end
end

#todo: pass in wrapper instead of redefining it each time based on input p
function local_newton_nestedtransition(F, x0, p; tol=1e-12, max_iter=30, under_relaxation=true, iters_direct=15)

    x = copy(x0)
    F_wrapper(x) = F(x, p)

    get_Hk_wrapper(states) = get_Hk(states[3], states[2], p.params; airfoil=true, turbulent=false, limit_Hk_max=false, limit_Hk_min=false)
    Hk_step = 0.0
    direct = true

    for i in 1:max_iter

        R = F_wrapper(x)
        # @infiltrate
        if norm(R) < tol # Convergence check
            return x, true
        end

        if direct
            J = ForwardDiff.jacobian(F_wrapper, x)
            dx = pinv(J) * -R
            # @infiltrate
        else
            J1 = ForwardDiff.jacobian(F_wrapper, x)
            J2 = ForwardDiff.gradient(get_Hk_wrapper, x)
            Jnew = vcat(J1, J2')
            Hk = get_Hk_wrapper(x)
            b = vcat(-R, Hk_step - Hk)

            dx = pinv(Jnew) * b
            # @infiltrate
        end

        #under relaxation
        if under_relaxation
            dm = max(abs(dx[2]/p.U1[1]), abs(dx[3]/p.U1[2]))

            dm = max(dm, abs(dx[3]/10))
            
            omega = dm > 0.3 ? 0.3/dm : 1.0
            # @infiltrate
            dx[2:3] .*= omega #TODO: this and the next line yield the same result for now - check on this later
            # dx .*= omega
        end

        dmax = 0.2 * p.dx * (1.1 - i/max_iter)

        if abs(dx[1]) > dmax
            # @infiltrate
            dx[1] = dx[1] * dmax / abs(dx[1])
        end

        xnew = x + dx

        ####### check for separation #######
        # Hk = get_Hk_wrapper(xnew)        
        # Hkmax = 3.8 #laminar

        # if direct && ((Hk > Hkmax) || i > iters_direct || any(xnew[1:2] .< 0)) 
        #     #don't take the update; prescribe Hk instead; skip if already separated
        #         @infiltrate
        #     direct = false
        #     Hk = get_Hk(p.U1[2], p.U1[1], p.params; airfoil=true, turbulent=false, limit_Hk_max=false, limit_Hk_min=false)
        #     Hkr = p.dx / p.U1[1]

        #     Hk_step = Hk + 0.03 * Hkr

        #     Hk_step = max(Hk_step, Hkmax)

        #     # if i > iters_direct && !debug
        #     #     x = copy(x0) #reinitialize states
        #     # end

        # else
        #     x = xnew #take the update
        # end
    end

    # check one more time
    R = F_wrapper(x)

    converged = false
    if norm(R) < tol # Convergence check
        converged = true
    end

    return x, converged
end

function local_newton_xt(F, x0, p; tol=1e-12, max_iter=20)

    x = copy(x0)
    F_wrapper(x) = F(x, p)

    for i in 1:max_iter

        R = F_wrapper(x)
        if norm(R) < tol # Convergence check
            return x, true
        end

        J = ForwardDiff.derivative(F_wrapper, x)

        dx = -R / J

        dmax = 0.2 * p.dx * (1.1 - i/max_iter)

        if abs(dx) > dmax
            dx = dx * dmax / abs(dx)
        end

        x = x + dx
    end

    # check one more time
    R = F_wrapper(x)

    if norm(R) < tol # Convergence check
        return x, true
    else
        return x, false
    end
end

function local_newton_implicit(F, x0; kwargs...)

    solution, converged = local_newton(F, x0; kwargs...)

    return solution
end

function local_newton_amp(F, x0; tol=1e-10, max_iter=50, debug=false, iters_direct=12)
    x = copy(x0)

    for i in 1:max_iter
        Fx = F(x)

        if norm(Fx) < tol # Convergence check
            return x, true
        end
            
        J = ForwardDiff.gradient(F, x)
        Δx = pinv(J) * -Fx  # Linear system solve

        #under relaxation
        dm = 0.5 * (1.01 - i/max_iter)
        omega = abs(Δx[1]) > dm ? dm / abs(Δx[1]) : 1.0
        Δx .*= omega

        x += Δx

        # Prevent extreme values
        if !F.turbulent
            if x[1] < 0
                x[1] = 0
            end
        end
    end

    # check one more time
    Fx = F(x)

    if norm(Fx) < tol # Convergence check
        return x, true
    else
        return x, false
    end
end

function global_newton(F, U0; tol=1e-10, max_iter=50, under_relaxation=true, verbose=true)
    F_wrapper(x) = F(x, nothing)

    U = copy(U0)
    num_nodes_total = Int(length(F.d))
    num_node_states = Int(num_nodes_total * 4) # excluding transition states

    # if F.operatingparameters.transition_method == 3
    #     idx_thetas = vcat(1:4:num_node_states, num_node_states+2, num_node_states+5)
    #     idx_deltas = vcat(2:4:num_node_states, num_node_states+3, num_node_states+6)
    # else
        idx_thetas = vcat(1:4:num_node_states)
        idx_deltas = vcat(2:4:num_node_states)
    # end
    idx_state3s = 3:4:num_node_states
    idx_ues = 4:4:num_node_states

    for i in 1:max_iter
        R = F_wrapper(U)
        R_norm = norm(R[1:3*length(idx_ues)]) #what mfoil uses, but doesn't seem right to exclude ue states
        # R_norm = norm(R)

        if verbose; println("\nNewton iteration $i: Rnorm = $(norm(R[1:3*length(idx_ues)]))"); end
        # println("\nNewton iteration $i: Rnorm = $(norm(R))")

        if R_norm < tol # Convergence check
            if verbose
                println("")
                println("Converged at step $(i)!")
                println("R_norm: ", R_norm)
                println("")
            end
            return U, true
        end

        # define states
        deltas = U[idx_deltas]
        thetas = U[idx_thetas]
        state3s = U[idx_state3s] 
        ues = U[idx_ues] 

        ctmax = maximum(state3s[F.turbulent_nodes])

        cj = ForwardDiff.JacobianConfig(F_wrapper, U, ForwardDiff.Chunk{4}())
        J = ForwardDiff.jacobian(F_wrapper, U, cj)  # Automatic Jacobian

        # solve system
        ΔU = pinv(J) * -R

        ##### Force ΔU to zero in certain cases to help stability #####
        for j in idx_state3s
            if U[j] == 0 && abs(ΔU[j]) < 1e-8
                ΔU[j] = 0.0
            end
        end

        Δdeltas = ΔU[idx_deltas]
        Δthetas = ΔU[idx_thetas]
        Δstate3s = ΔU[idx_state3s] 
        Δues = ΔU[idx_ues] 
        
        omega = 1.0

        # prevent big decrease in theta
        fmin = minimum(Δthetas ./ thetas)
        om = fmin < -0.5 ? abs(0.5 / fmin) : 1.0
        omega = om < omega ? om : omega

        # prevent big decrease in deltastar
        fmin = minimum(Δdeltas ./ deltas)
        om = fmin < -0.5 ? abs(0.5 / fmin) : 1.0
        omega = om < omega ? om : omega

        # limit negative n/ctau #TODO: is there a simpler way to code this out?
        for j in 1:length(state3s)
            if !F.turbulent_nodes[j] && state3s[j] < 0.2
                nothing
            elseif F.turbulent_nodes[j] && state3s[j] < 0.1 * ctmax
                nothing
            elseif state3s[j] == 0.0 || Δstate3s[j] == 0.0
                nothing
            elseif state3s[j] + Δstate3s[j] < 0.0
                om = 0.8 * abs(state3s[j] / Δstate3s[j])
                omega = om < omega ? om : omega
            end
        end

        # prevent large changes in n
        dnmax = maximum(abs.(Δstate3s[.!F.turbulent_nodes]))
        om = dnmax > 0.0 ? abs(2 / dnmax) : 1.0
        omega = om < omega ? om : omega

        # prevent large changes in ctau
        dctmax = maximum(abs.(Δstate3s[F.turbulent_nodes]))
        om = dctmax > 0.0 ? abs(0.5 / dctmax) : 1.0
        omega = om < omega ? om : omega

        # prevent large changes in ue
        fmax = maximum(abs.(Δues) / F.operatingparameters.Vinf)
        om = fmax > 0.0 ? abs(0.2 / fmax) : 1.0
        omega = om < omega ? om : omega

        # take the update
        U += ΔU * omega

        # fix bad Hk after update
        for j in 1:length(thetas)
            if j < F.idx_w[1]
                Hkmin = 1.02
            else
                Hkmin = 1.00005
            end
            Hk = get_Hk(U[idx_deltas[j]], U[idx_thetas[j]], U[idx_ues[j]], F.operatingparameters; limit_Hk_max=false, limit_Hk_min=false)

            if Hk < Hkmin
                U[idx_deltas[j]] += 2 * (Hkmin - Hk) * U[idx_thetas[j]]
            end
        end

        # fix negative ctau after update
        for j in 1:length(state3s)
            if F.turbulent_nodes[j] && U[idx_state3s[j]] < 0
                U[idx_state3s[j]] = 0.1 * ctmax
            end
        end

        ################ Move Stagnation Point ################
        _, _, _, _, ues_new = stagnation_point_move!(F, U[idx_ues])
        if any(ues_new .< 0)
            @infiltrate
        end
        U[idx_ues] .= ues_new

        if F.operatingparameters.transition_method == 3
            U_transition = U[end-5:end]
            U_temp = reshape(U[1:end-6], 4, :)'

            ################ Re-calculate Transition ################
            resolve_transition_location!(F, U_temp; upper=false)
            resolve_transition_location!(F, U_temp; upper=true)

            U_temp = reshape(U_temp', :, 1) #todo: probably a better way to do this...
            U = vcat(U_temp, U_transition)
        else
            U = reshape(U, 4, :)'

            ################ Re-calculate Transition ################
            resolve_transition_location!(F, U; upper=false)
            try
                resolve_transition_location!(F, U; upper=true)
            catch
                @infiltrate
            end

            U = reshape(U', :, 1) #todo: probably a better way to do this...
        end
    end

    # check residual one more time
    R = F_wrapper(U)
    # R_norm = norm(R)
    R_norm = norm(R[1:3*length(idx_ues)])

    if R_norm < tol # Convergence check
        if F.operatingparameters.verbose
            println("Converged at last step!")
            println("R_norm: ", R_norm)
        end
        return U, true
    else
        if F.operatingparameters.verbose
            println("Failed to converge:")
            println("R_norm: ", R_norm)
        end
        return U, false
    end
end

function local_levenberg_marquardt(F, x0; tol=1e-10, max_iter=100, λ_init=1e-3, under_relaxation=false)
    x = copy(x0)
    λ = λ_init

    for i in 1:max_iter

        res = F(x) #todo: rearrange so it only calls once per iteration (move this out in front of loop, redefine based on res_new)

        # Convergence check
        if norm(res) < tol
            return x, true
        end

        J = ForwardDiff.jacobian(F, x)
        A = J' * J + λ * I
        g = J' * res

        Δx = pinv(A) * -g  # Linear system solve

        #under relaxation
        if under_relaxation
            dm = max(abs(Δx[1]/x[1]), abs(Δx[2]/x[2]))
            omega = dm > 0.3 ? 0.3/dm : 1.0
            Δx .*= omega
        end

        x_new = x + Δx

        # Prevent extreme values
        if F.turbulent
            if x_new[3] > 0.3
                x_new[3] = 0.3
            elseif x_new[3] < 1e-7
                x_new[3] = 1e-7
            end
        else
            if x_new[3] < 0
                x_new[3] = 0
            end
        end

        res_new = F(x_new)

        # Accept or reject step
        if norm(res_new) < norm(res)
            x = x_new
            λ /= 10
            #todo: res = res_new
        else
            λ *= 10
        end
    end

    # Check convergence one more time
    res = F(x)

    if norm(res) < tol
        return x, true
    else
        return x, false
    end
end

function global_levenberg_marquardt(F, U0; tol=1e-10, max_iter=100, λ_init=1e-3)
    F_wrapper(x) = F(x, nothing)
    
    U = copy(U0)
    λ = λ_init

    num_nodes_total = Int(length(F.d))
    idx_deltas = 1:num_nodes_total
    idx_thetas = num_nodes_total+1:num_nodes_total*2
    idx_state3s = num_nodes_total*2+1:num_nodes_total*3
    idx_ues = num_nodes_total*3+1:num_nodes_total*4

    for i in 1:max_iter

        Fx = F_wrapper(U) 

        # Convergence check
        if norm(Fx) < tol
            return U, true
        end

        # define states
        deltas = U[idx_deltas]
        thetas = U[idx_thetas]
        state3s = U[idx_state3s] 
        ues = U[idx_ues] 

        ctmax = maximum(state3s[F.turbulent_nodes])

        J = ForwardDiff.jacobian(F_wrapper, U)
        A = J' * J + λ * I
        g = J' * Fx

        ΔU = pinv(A) * -g  # Linear system solve

        Δdeltas = ΔU[1:num_nodes_total]
        Δthetas = ΔU[idx_thetas]
        Δstate3s = ΔU[idx_state3s] 
        Δues = ΔU[idx_ues] 
        
        omega = 1.0

        # First prevent big decrease in deltastar
        fmin = minimum(Δdeltas ./ deltas)
        om = fmin < -0.5 ? abs(0.5 / fmin) : 1.0
        omega = om < omega ? om : omega

        # prevent big decrease in theta
        fmin = minimum(Δthetas ./ thetas)
        om = fmin < -0.5 ? abs(0.5 / fmin) : 1.0
        omega = om < omega ? om : omega

        # limit negative n/ctau #TODO: is there a simpler way to code this out?
        for i in 1:length(state3s)
            if !F.turbulent_nodes[i] && state3s[i] < 0.2
                nothing
            elseif F.turbulent_nodes[i] && state3s[i] < 0.1 * ctmax
                nothing
            elseif state3s[i] == 0.0 || Δstate3s[i] == 0.0
                nothing
            elseif state3s[i] + Δstate3s[i] < 0.0
                om = 0.8 * abs(state3s[i] / Δstate3s[i])
                omega = om < omega ? om : omega
            end
        end

        # prevent large changes in n
        dnmax = maximum(abs.(Δstate3s[.!F.turbulent_nodes]))
        om = dnmax > 0.0 ? abs(2 / dnmax) : 1.0
        omega = om < omega ? om : omega

        # prevent large changes in ctau
        dctmax = maximum(abs.(Δstate3s[F.turbulent_nodes]))
        om = dctmax > 0.0 ? abs(0.5 / dctmax) : 1.0
        omega = om < omega ? om : omega

        # prevent large changes in ue
        fmax = maximum(abs.(Δues) / F.operatingparameters.Vinf)
        om = fmax > 0.0 ? abs(0.2 / fmax) : 1.0
        omega = om < omega ? om : omega

        U_new = U + ΔU * omega

        # fix bad Hk after update
        for i in 1:length(thetas)
            if i < F.idx_w[1]
                Hkmin = 1.02
            else
                Hkmin = 1.00005
            end
            Hk = get_Hk(U[idx_deltas[i]], U[idx_thetas[i]], U[idx_ues[i]], F.operatingparameters; limit_Hk_max=false, limit_Hk_min=false)

            if Hk < Hkmin
                U[idx_deltas[i]] += 2 * (Hkmin - Hk) * U[idx_thetas[i]]
            end
        end

        # fix negative ctau after update
        for i in 1:length(state3s)
            if F.turbulent_nodes[i] && U[idx_state3s[i]] < 0
                U[idx_state3s[i]] = 0.1 * ctmax
            end
        end

        # fix negative n after update
        for i in 1:length(state3s)
            if !F.turbulent_nodes[i] && U[idx_state3s[i]] < 0
                U[idx_state3s[i]] = 0.0
            end
        end

        Fx_new = F_wrapper(U_new)

        # Accept or reject step
        if norm(Fx_new) < norm(Fx)
            U = U_new
            λ /= 10
        else
            λ *= 10
        end
    end

     # Check convergence one more time
     Fx = F_wrapper(U)

     if norm(Fx) < tol
         return U, true
     else
        println("")
        println("Final Global Residual norm: $(norm(Fx))")
        println("")
         return U, false
     end
end
