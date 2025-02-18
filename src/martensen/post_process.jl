function post_process(
    method::Martensen, problem, panels::TP, mesh, solution; npanels=80, debug=false
) where {TP<:Panel}
    return post_process(
        method::Martensen, problem, [panels], mesh, solution; npanels=80, debug=false
    )
end

function post_process(
    ::Martensen, problem, panels, mesh, solution; npanels=80, debug=false
)

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    nbodies = mesh.nbodies
    flow_angle = problem.flow_angle
    naoa = length(flow_angle)

    # - Initialize Outputs - #
    TF = eltype(mesh.x)
    # Surface Velocities
    v_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    # Surface Pressures
    p_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    xsmooth = zeros(TF, nbodies, 2 * npanels - 1, naoa)

    # vortex strengths
    gamma0 = solution.x[:, 1]
    gamma90 = solution.x[:, 2]
    gamma_u = 0.0
    gamma_v = 0.0

    #compute unit solutions gamma_u and gamma_v
    for i = 1:nbodies
        gamma_u = gamma_u + x[i, 1]*panels.panel_length[i]
        gamma_v = gamma_v + x[i, 2]*panels.panel_length[i]
    end

    # Flow angles
    flow_angle = problem.flow_angle
    
    #compute solution parameters
    pitch = 1.0 #not sure were pitch is going to be input yet
    k1 = (1 - gamma_v / (2 * pitch)) / (1 + gamma_v / (2 * pitch))
    k2 = gamma_u / (pitch*(1 + gamma_v / (2*pitch)))
    beta2 = atan(k1*tan(flow_angle) - k2)
    betainf = atan(0.5*(sin(flow_angle) / cos(flow_angle) + sin(beta2) / cos(beta2)))
    Winf = W1*cos(flow_angle) / cos(betainf)
    Uinf = Winf*cos(betainf)
    Vinf = Winf*sin(betainf)

    for m in 1:nbodies
        for a in 1:naoa
            # - Extract surface velocity - #
            vti = solution.x[1:idx[end][end]]
            # vti = [
            #     gamma0[i] * cosd(flow_angle[a]) + gamma90[i] * sind(flow_angle[a]) for
            #     i in idx[m]
            # ]

            # - Calculate surface pressure - #
            cpi = 1.0 .- (vti) .^ 2

            ### --- Smooth Distributions --- ###
            #smooth_distributions functions are found in utils.jl
            #smooth_distributions functions are found in utils.jl
            v_surf[m, :, a], xsmooth[m, :, a] = smooth_distributions(
                Constant(), panels[m].panel_center, vti, npanels
            )
            p_surf[m, :, a], _ = smooth_distributions(
                Constant(), panels[m].panel_center, cpi, npanels
            )
        end
    end

    return PeriodicPolar(v_surf, p_surf, xsmooth)
end
