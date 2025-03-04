function post_process(method::HessSmith, panel_geometry, system_geometry, strengths, flow_angles)
    return post_process(
        method, [panel_geometry], system_geometry, strengths, flow_angles
    )
end

"""

Calculate tangential velocity and surface pressures

"""
function post_process(
    method::HessSmith, panel_geometry::AbstractVector, system_geometry, strengths, flow_angles
)

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    nbodies = system_geometry.nbodies
    naoa = length(flow_angles)

    # - Initialize Outputs - #
    TF = eltype(system_geometry.r_influence)

    tangential_velocities = [
        zeros(idx[m][end] - idx[m][1] + 1, length(flow_angles)) for m in 1:nbodies
    ]
    surface_pressures = [
        zeros(idx[m][end] - idx[m][1] + 1, length(flow_angles)) for m in 1:nbodies
    ]
    summed_strengths = zeros(maximum(idx[m][end] for m in 1:nbodies), naoa)

    for a in 1:naoa
        strengths[:, 1] *= cos(flow_angles[a])
        strengths[:, 2] *= sin(flow_angles[a])
        summed_strengths[:, a] = sum(strengths; dims=2) * 2 * pi * method.V_inf
    end

    for m in 1:nbodies  # Loop through number of bodies
        for a in 1:naoa # Loop through the different angles of attack         
            for i in idx[m][1]:idx[m][end] # Loop through the panels
                set1 = 0.0
                set2 = 0.0
                for j in idx[m][1]:idx[m][end]-1   
                    set1 += summed_strengths[j, a] * (system_geometry.beta[i, j] * system_geometry.sine_angle_panels[i, j] - log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.cos_angle_panels[i, j])
                    set2 += system_geometry.beta[i, j] * system_geometry.cos_angle_panels[i, j] + log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.sine_angle_panels[i, j]
                end
                tangential_velocities[m][i - idx[m][1] + 1, a] = method.V_inf * (panel_geometry.cosine_vector[i] * cos(flow_angles[a]) + panel_geometry.sine_vector[i] * sin(flow_angles[a])) + (set1 / (2 * π)) + (summed_strengths[end, a] / (2 * π)) * set2
                surface_pressures[m][i - idx[m][1] + 1, a] = 1.0 - tangential_velocities[m][i - idx[m][1] + 1, a] ^ 2
            end
        end
    end

    return (; tangential_velocities, surface_pressures)
end

# Loop bodies, angles of attack, and panels