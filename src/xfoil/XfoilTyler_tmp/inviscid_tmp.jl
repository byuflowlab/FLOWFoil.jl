
#### inviscid functions that maybe shadow what the original code does. Necessary for B, C, D matrices ####
function inviscid_velocities_airfoil_and_wake(sinalpha, cosalpha, γref, ue_wake_ref)

    #TODO: import sign directions for all nodes (when initializing won't need it, but may later?)
    cs = [cosalpha, sinalpha]  
    
    Vairfoil = abs.(γref * cs)
    Vwake = ue_wake_ref * cs

    Vwake[1] = Vairfoil[end]

    return vcat(Vairfoil, Vwake)
end

function inviscid_velocities_airfoil(sinalpha, cosalpha, γref)

    #TODO: import sign directions for all nodes (when initializing won't need it, but may later?)
    Vairfoil = γref * [cosalpha, sinalpha] 

    return Vairfoil
end

function get_invisicd_velocity_from_gammas(panels_airfoil, gammas, xbar, panel_TE, operatingparameters; cosalpha=nothing, sinalpha=nothing)

    xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[end].x2)
    ys = vcat([p.y1 for p in panels_airfoil], panels_airfoil[end].y2)
    lengths = [p.length for p in panels_airfoil]
    cthetas = [p.ctheta for p in panels_airfoil]
    sthetas = [p.stheta for p in panels_airfoil]

    Vinf = operatingparameters.Vinf
    if isnothing(cosalpha)
        cosalpha = operatingparameters.cosalpha[1] #todo: edit for multiple alphas
        sinalpha = operatingparameters.sinalpha[1]
    end

    V_influence = zeros(2)
    # loop over airfoil panels for linear vortex influence
    for j in 1:length(panels_airfoil)

        xstar, ystar, _, _, 
            logr1, logr2, beta = get_relative_geometry_panel_point(xbar[1], xbar[2], 
                                                                    xs[j], xs[j+1], 
                                                                    ys[j], ys[j+1], 
                                                                    cthetas[j], sthetas[j], 
                                                                    lengths[j])

        a, b = get_linear_vortex_velocity_influence(xstar, ystar, logr1, logr2, 
                                                    beta, lengths[j], xs[j], ys[j], 
                                                    xs[j+1], ys[j+1])
        V_influence .+= a * gammas[j] + b * gammas[j+1]
    end

    # TE source influence
    xstar_TE, ystar_TE, _, _, 
            logr1_TE, logr2_TE, beta_TE = get_relative_geometry_panel_point(xbar[1], xbar[2], 
                                                                    panel_TE.x1, panel_TE.x2, 
                                                                    panel_TE.y1, panel_TE.y2,
                                                                    panel_TE.ctheta, panel_TE.stheta, 
                                                                    panel_TE.length)

    a_source = get_constant_source_velocity_influence(xstar_TE, ystar_TE, panel_TE.length, 
                                                logr1_TE, logr2_TE, beta_TE, 
                                                panel_TE.x1, panel_TE.y1, 
                                                panel_TE.x2, panel_TE.y2)

    f_source = a_source * panel_TE.tcp / 2
    V_influence .+= f_source * gammas[end] - f_source * gammas[1]

    # TE vortex influence
    a_vortex, b_vortex = get_linear_vortex_velocity_influence(xstar_TE, ystar_TE, logr1_TE, logr2_TE, 
                                                    beta_TE, panel_TE.length, panel_TE.x1, panel_TE.y1, 
                                                    panel_TE.x2, panel_TE.y2)
    f_vortex = (a_vortex + b_vortex) * panel_TE.tdp / 2
    V_influence .+= f_vortex * gammas[1] - f_vortex * gammas[end]

    # freestream influence
    V_influence .+= Vinf * [cosalpha, sinalpha]
    return V_influence
end

function get_linear_vortex_velocity_influence(xstar, ystar, logr1, logr2, beta, lj, x1, y1, x2, y2)

    t = [x2-x1; y2-y1]; t = t/norm(t);
    n = [-t[2]; t[1]]

    return get_linear_vortex_velocity_influence(xstar, ystar, logr1, logr2, beta, lj, t, n)
end

function get_linear_vortex_velocity_influence(xstar, ystar, logr1, logr2, beta, lj, t, n)

    logr1mr2 = logr1 - logr2

    Q1 = beta / (2pi)
    Q2 = (ystar * logr1mr2 - xstar * beta) / (2pi * lj)

    u1 = Q1 + Q2
    u2 = -Q2

    Q3 = -logr1mr2 / (2pi)
    Q4 = (xstar * logr1mr2 - lj + ystar * beta) / (2pi * lj)
    
    v1 = Q3 + Q4
    v2 = -Q4

    a = [u1 * t[1] + v1 * n[1]; u1 * t[2] + v1 * n[2]]
    b = [u2 * t[1] + v2 * n[1]; u2 * t[2] + v2 * n[2]]

    return a, b #a is multiplied by gamma, b by gamma+1
end


function get_linear_source_velocity_influence(xstar, ystar, logr1, logr2, beta, lj, x1, y1, x2, y2; vdir=nothing)

    t = [x2-x1; y2-y1]; t = t/norm(t);
    n = [-t[2]; t[1]]

    return get_linear_source_velocity_influence(xstar, ystar, logr1, logr2, beta, lj, t, n; vdir)
end

function get_linear_source_velocity_influence(xstar, ystar, logr1, logr2, beta, lj, t, n; vdir=nothing)

    logr1mr2 = logr1 - logr2

    Q1 = logr1mr2 / (2pi)
    Q2 = (xstar * logr1mr2 - lj + ystar * beta) / (2pi * lj)
    
    u1 = Q1 - Q2
    u2 = Q2
    
    Q3 = beta / (2pi)
    Q4 = (xstar * beta - ystar * logr1mr2) / (2pi * lj)
    
    v1 = Q3 - Q4
    v2 = Q4

    a = [u1 * t[1] + v1 * n[1]; u1 * t[2] + v1 * n[2]]
    b = [u2 * t[1] + v2 * n[1]; u2 * t[2] + v2 * n[2]]

    if !isnothing(vdir)
        a = a' * vdir
        b = b' * vdir
    end
    return a, b #a is multiplied by sigma, b by sigma+1
end

function get_constant_source_velocity_influence(xstar, ystar, lj, logr1, logr2, β, x1, y1, x2, y2; vdir=nothing)

    t = [x2-x1; y2-y1]; t = t/norm(t);
    n = [-t[2]; t[1]]

    return get_constant_source_velocity_influence(xstar, ystar, lj, logr1, logr2, β, t, n; vdir)
end

function get_constant_source_velocity_influence(xstar, ystar, lj, logr1, logr2, β, t, n; vdir=nothing)

    theta1 = atan(ystar, xstar)
    theta2 = atan(ystar, xstar-lj)

    u = 1/(2π) * (logr1 - logr2)
    v = 1/(2π) * (theta2 - theta1)   

    a = [u*t[1] + v*n[1]; u*t[2] + v*n[2]]
    if !isnothing(vdir)
        a = a' * vdir
    end

    return a
end

function get_constant_vortex_streamfunction_influence(I1, t_TE, p_TE)
    return I1/2 * dot(t_TE, p_TE) 
end

function get_constant_vortex_streamfunction_influence(xstar, ystar, logr, logrp1, lj, β, t_TE, p_TE)

    I1 = get_I1(xstar, ystar, logr, logrp1, lj, β)

    return I1/2 * abs(dot(t_TE, p_TE))
end


function get_constant_source_streamfunction_influence(xstar, ystar, logr, logrp1, lj, β)

    theta1 = atan(ystar, xstar)
    theta2 = atan(ystar, xstar-lj)

    if logr == 0.0
        # if iszero(logr)
        theta1 = pi
        theta2 = pi
        β = 0
    elseif logrp1 == 0.0
        # elseif iszero(logrp1)
        theta1 = 0
        theta2 = 0
        β = 0
    end

    I3 = get_I3(xstar, ystar, logr, logrp1, lj, β, theta2)

    #shift source influence to ensure branch cut orientation is behaving well
    if theta1 + theta2 > pi #todo redo in terms of theta and beta?
        #todo: make sure this is ok, might need to add a abs() or say < -pi
        I3 -= 0.25 * lj
    else
        I3 += 0.75 * lj
    end

    return I3
end

function get_linear_source_streamfunction_influence(xstar, ystar, logr, logrp1, lj, β, r, rp1)

    I3 = get_I3(xstar, ystar, logr, logrp1, lj, β; linear=true)
    I5 = get_I5(xstar, ystar, lj, β, r, rp1)

    a = I3 - (xstar*I3 + I5) / lj # to be multiplied by sigma_k
    b = (xstar*I3 + I5) / lj

    return a, b
end

function get_I1(xstar, ystar, logr, logrp1, lj, β)

    I1 = 1/(2pi) * (xstar * (logr - logrp1) + lj*logrp1 + ystar*β - lj)

    return I1
end

function get_I3(xstar, ystar, logr, logrp1, lj, β, theta2=atan(ystar, xstar-lj); linear=false, debug=false) #for constant source term

    # theta = atan(xstar - lj, ystar)
    # theta = -atan(ystar, xstar - lj)
    # theta = atan(xstar - lj, ystar) # if I do this route, I need to add 2pi I think if it's negative to make a branch cut
    # theta2 = atan(ystar, xstar-lj)

    if linear && theta2 < 0 #! see if this is causing issues
        theta2 += 2pi
    end

    I3 = 1/(2pi) * (ystar*(logr - logrp1) - xstar * β + lj * theta2)
    # I3 = 1/(2pi) * (ystar*(logr - logrp1) - xstar * β + lj * (π/2 - theta))
    # I3 = 1/(2pi) * (ystar*(logrp1 - logr) + xstar * β + lj * theta)

    return I3
end


function get_I5(xstar, ystar, lj, β, r, rp1)

    # theta = atan(xstar - lj, ystar)
    theta2 = atan(ystar, xstar-lj)
    theta1 = atan(ystar, xstar)

    #branch cuts
    if theta1 < 0; theta1 += 2pi; end
    if theta2 < 0; theta2 += 2pi; end

    I5 = 1 / (4pi) * (rp1^2*theta2 - r^2*theta1 - ystar*lj) #! MFOIL math
    #TODO: match these in the theory
    # I5 = 1 / (4pi) * (ystar * lj + β * r^2 + theta * (r^2 - rp1^2)) #! My math - incorrect

    return I5
end
#####################################################################