function post_process(
    method::Martensen, panel_geometry, system_geometry, strengths, flow_angles
) where {TP<:Panel}
    return post_process(
        method::Martensen, [panel_geometry], system_geometry, strengths, flow_angles
    )
end

function post_process(
    method::Martensen,
    panel_geometry::AbstractVector,
    system_geometry,
    strengths,
    flow_angles,
)

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    nbodies = system_geometry.nbodies
    naoa = length(flow_angles)

    # - Initialize Outputs - #
    TF = eltype(system_geometry.x)
    # Surface Velocities
    tangential_velocities = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    surface_pressures = zeros(TF, nbodies, 2 * npanels - 1, naoa)

    # vortex strengths
    gamma0 = strengths[:, 1]
    gamma90 = strengths[:, 2]

    #TODO:TODO:TODO:TODO:TODO:TODO:TODO:TODO>>>>>
    #WHAT IS THIS???
    gamma_u = 0.0
    gamma_v = 0.0
    #compute unit solutions gamma_u and gamma_v
    for i in 1:nbodies
        #WHAT'S GOING ON HERE???
        gamma_u = gamma_u + x[i, 1] * panel_geometry.panel_length[i]
        gamma_v = gamma_v + x[i, 2] * panel_geometry.panel_length[i]
    end
    #<<<<<TODO:TODO:TODO:TODO:TODO:TODO:TODO:TODO

    #compute strengths parameters
    pitch = method.pitch
    W1 = 1.0 #freestream velocity
    #TODO:TODO:TODO:TODO:TODO:TODO:TODO:TODO>>>>>
    #should gamma_u and gamma_v actually be gamma0 and gamma90, respectively?
    #<<<<<TODO:TODO:TODO:TODO:TODO:TODO:TODO:TODO
    k1 = (1 - gamma_v / (2 * pitch)) / (1 + gamma_v / (2 * pitch))
    k2 = gamma_u / (pitch * (1 + gamma_v / (2 * pitch)))

    for m in 1:nbodies
        for a in 1:naoa

            # - Calculate Velocity Components - #
            beta2 = atan(k1 * tan(flow_angles[a]) - k2)
            betainf = atan(
                0.5 * (sin(flow_angles[a]) / cos(flow_angles[a]) + sin(beta2) / cos(beta2))
            )
            Winf = W1 * cos(flow_angles[a]) / cos(betainf)
            Uinf = Winf * cos(betainf)
            Vinf = Winf * sin(betainf)

            # - Calculate surface velocity - #
            tangential_velocities[m, :, a] = [
                (Uinf * gamma0[i] + Vinf * gamma90[i]) / W1 for i in idx[m]
            ]

            # - Calculate surface pressure - #
            surface_pressures[m, :, a] = 1.0 .- (tangential_velocities) .^ 2
        end
    end

    return (; tangential_velocities, surface_pressures)
end
