function post_process(
    ap::PeriodicProblem, problem, panels::TP, mesh, solution; npanels=80, debug=false
) where {TP<:Panel}
    return post_process(
        ap::PeriodicProblem, problem, [panels], mesh, solution; npanels=80, debug=false
    )
end

function post_process(
    ::PeriodicProblem, problem, panels, mesh, solution; npanels=80, debug=false
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

    # Flow angles
    flow_angle = problem.flow_angle

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
