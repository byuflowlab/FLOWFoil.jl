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
    # Convert to radians
    flow_angles = (pi / 180.0) * flow_angles

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    nbodies = system_geometry.nbodies
    naoa = length(flow_angles)

    # - Initialize Outputs - #
    TF = eltype(system_geometry.r_x)

    vs = [
        zeros(idx[m][end] - idx[m][1] + 1, naoa) for m in 1:nbodies
    ]
    cp = [
        zeros(idx[m][end] - idx[m][1] + 1, naoa) for m in 1:nbodies
    ]

    cl = zeros(naoa, nbodies)
    cd = zeros(naoa, nbodies)
    cm = zeros(naoa, nbodies)

    for m in 1:nbodies
        # vortex strengths per unit length
        gamma0 = [strengths[idx[m][1:(end - 1)]  .- nbodies*(m-1), 1]; -strengths[idx[m][1] .- nbodies*(m-1) + 1, 1]]
        gamma90 = [strengths[idx[m][1:(end - 1)] .- nbodies*(m-1), 2]; -strengths[idx[m][1] .- nbodies*(m-1) + 1, 2]]

        # total circulation
        gamma_u = dot(panel_geometry[m].panel_length, gamma0)
        gamma_v = dot(panel_geometry[m].panel_length, gamma90)

        #compute strengths parameters

        Vinf = 1.0 #freestream velocity

        if method.cascade
            k1 =
                (1.0 - gamma_v / 2.0 / system_geometry.pitch) /
                (1.0 + gamma_v / 2.0 / system_geometry.pitch)
            k2 =
                gamma_u / system_geometry.pitch /
                (1.0 + gamma_v / 2.0 / system_geometry.pitch)
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
            if m == 1
                vs[m][:, a] = [
                    (w_x * gamma0[i] + w_y * gamma90[i]) / Vinf for i in idx[m]
                ]
            else
                vs[m][:, a] = [
                    (w_x * gamma0[i] + w_y * gamma90[i]) / Vinf for i in idx[m] .- idx[m-1][end]
                ]
            end

            # - Calculate surface pressure - #
            cp[m][:, a] = 1.0 .- vs[m][:, a] .^ 2

            #compute cascade lift if cascade model is true, else use planar lift
            if method.cascade
                cl[a,m] =
                    2.0 *
                    system_geometry.pitch *
                    (tan(flow_angles[a]) - tan(beta2)) *
                    cos(betainf) #Important: This equation assumes that the chord length is 1.0
            else                
                cl[a,m] =
                    2.0 * (gamma_u * w_x + gamma_v * w_y) /
                    (W * calculate_chord(panel_geometry))
                #=
                lift_coefficients[m][a] = 2*(cos(flow_angles[a])*sum(gamma0 .* panel_geometry[m].panel_length)
                + sin(flow_angles[a])*sum(gamma90 .* panel_geometry[m].panel_length)) #Important: This equation assumes that the chord length is 1.0
                =#
            end
        end
    end
    if nbodies == 1
        #if it is a single body, this reduces the need to use the body index
        vs_new = zeros(idx[1][end]-idx[1][1]+1, naoa)
        cp_new = zeros(idx[1][end]-idx[1][1]+1, naoa)

        vs_new[:,:] = vs[1][:,:]
        cp_new[:,:] = cp[1][:,:]
        vs = vs_new
        cp = cp_new
        cl = cl[:]
        cd = cd[:]
        cm = cm[:]
    end

    return (; vs, cp, cl, cd, cm)
end
