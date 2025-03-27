function post_process(
    method::Martensen, panel_geometry, system_geometry, strengths, flow_angles
)
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
    TF = eltype(system_geometry.r_x)

    tangential_velocities = [
        zeros(idx[m][end] - idx[m][1] + 1, length(flow_angles)) for m in 1:nbodies
    ]
    surface_pressures = [
        zeros(idx[m][end] - idx[m][1] + 1, length(flow_angles)) for m in 1:nbodies
    ]
    lift_coefficients = [
        zeros(1, length(flow_angles)) for m in 1:nbodies
    ]

    for m in 1:nbodies

        # vortex strengths per unit length
        # gamma0 = strengths[:, 1]
        gamma0 = [strengths[idx[m][1:(end - 1)], 1]; -strengths[idx[m][1], 1]]
        # gamma90 = strengths[:, 2]
        gamma90 = [strengths[idx[m][1:(end - 1)], 2]; -strengths[idx[m][1], 2]]

        # total circulation
        gamma_u = dot(panel_geometry[m].panel_length, gamma0)
        gamma_v = dot(panel_geometry[m].panel_length, gamma90)

        #compute strengths parameters
        Vinf = 1.0 #freestream velocity
        if method.cascade
            k1 = (1.0 - gamma_v / 2.0 / method.pitch) / (1.0 + gamma_v / 2.0 / method.pitch)
            k2 = gamma_u / method.pitch / (1.0 + gamma_v / 2.0 / method.pitch)
        else
            k1 = 1.0
            k2 = 0.0
        end

        for a in 1:naoa

            # - Calculate Velocity Components - #
            beta2 = atan(k1 * tan(flow_angles[a]) - k2)
            betainf = atan(0.5 * (tan(flow_angles[a]) + tan(beta2)))
            W = Vinf * cos(flow_angles[a]) / cos(betainf)
            w_x = W * cos(betainf)
            w_y = W * sin(betainf)

            # - Calculate surface velocity - #
            tangential_velocities[m][:, a] = [
                (w_x * gamma0[i] + w_y * gamma90[i]) / Vinf for i in idx[m]
            ]

            # - Calculate surface pressure - #
            surface_pressures[m][:, a] = 1.0 .- tangential_velocities[m][:, a] .^ 2

            #compute cascade lift if cascade model is true, else use planar lift
            if method.cascade
                cl_cascade = 2*method.pitch*(tan(flow_angles[a]) - tan(beta2))*cos(betainf) #Important: This equation assumes that the chord length is 1.0
                cl_planar = 2*(cos(flow_angles[a])*sum(gamma0 .* panel_geometry.panel_length)
                    + sin(flow_angles[a])*sum(gamma90 .* panel_geometry.panel_length))
                lift_coefficients[m][a] = FLOWMath.sigmoid_blend(
                    cl_cascade,
                    cl_planar,
                    method.solidity,
                    method.transition_value,
                    method.transition_hardness
                )
            else
                lift_coefficients[m][a] = 2*(cos(flow_angles[a])*sum(gamma0 .* panel_geometry.panel_length)
                + sin(flow_angles[a])*sum(gamma90 .* panel_geometry.panel_length)) #Important: This equation assumes that the chord length is 1.0
            end
        end
    end

    return (; tangential_velocities, surface_pressures, lift_coefficients)
end
