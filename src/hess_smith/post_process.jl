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

    # tangential_velocities = [
    #     zeros(idx[m][end] - idx[m][1] + 1, length(flow_angles)) for m in 1:nbodies
    # ]
    # pressure_coefficient = [
    #     zeros(idx[m][end] - idx[m][1] + 1, length(flow_angles)) for m in 1:nbodies
    # ]
    # summed_strengths = zeros((maximum(idx[m][end] for m in 1:nbodies) + 1), naoa)

    # # for a in 1:naoa
    # #     strengths[:, 1] *= cos(flow_angles[a])
    # #     strengths[:, 2] *= sin(flow_angles[a])
    # #     summed_strengths[:, a] = sum(strengths; dims=2) * 2 * pi * method.V_inf
    # # end
    # for a in 1:naoa
    #     strengths_copy = copy(strengths)
    #     strengths_copy[:, 1] *= cosd(flow_angles[a])
    #     strengths_copy[:, 2] *= sind(flow_angles[a])
    #     summed_strengths[:, a] = sum(strengths_copy; dims=2) * 2 * pi * method.V_inf
    # end

    # for m in 1:nbodies  # Loop through number of bodies
    #     for a in 1:naoa # Loop through the different angles of attack         
    #         for i in idx[m][1]:idx[m][end] # Loop through the panels
    #             set1 = 0.0
    #             set2 = 0.0
    #             for j in idx[m][1]:idx[m][end]-1   
    #                 set1 += summed_strengths[j, a] * (system_geometry.beta[i, j] * system_geometry.sine_angle_panels[i, j] - log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.cos_angle_panels[i, j])
    #                 set2 += system_geometry.beta[i, j] * system_geometry.cos_angle_panels[i, j] + log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.sine_angle_panels[i, j]
    #             end
    #             tangential_velocities[m][i - idx[m][1] + 1, a] = method.V_inf * (panel_geometry[m].cosine_vector[i] * cosd(flow_angles[a]) + panel_geometry[m].sine_vector[i] * sind(flow_angles[a])) + (set1 / (2 * π)) + (summed_strengths[end, a] / (2 * π)) * set2
    #             # surface_pressures[m][i - idx[m][1] + 1, a] = 1.0 - tangential_velocities[m][i - idx[m][1] + 1, a] ^ 2
    #             pressure_coefficient[m][i-idx[m][1] + 1, a] = 1.0 - (tangential_velocities[m][i - idx[m][1] + 1, a] / method.V_inf)^2
    #         end
    #     end
    # end

    tangential_velocities = [
        zeros(length(flow_angles), idx[m][end] - idx[m][1] + 1) for m in 1:nbodies
    ]
    pressure_coefficient = [
        zeros(length(flow_angles), idx[m][end] - idx[m][1] + 1) for m in 1:nbodies
    ]
    
    summed_strengths = zeros(naoa, (maximum(idx[m][end] for m in 1:nbodies) + 1))

    for m in 1:nbodies  # Loop through number of bodies
        for a in 1:naoa  # Loop through the different angles of attack         
            strengths_copy = copy(strengths)
            strengths_copy[:, 1] *= cosd(flow_angles[a])
            strengths_copy[:, 2] *= sind(flow_angles[a])
            summed_strengths[a, :] = sum(strengths_copy; dims=2) * 2 * pi * method.V_inf
            for i in idx[m][1]:idx[m][end]  # Loop through the panels
                set1 = 0.0
                set2 = 0.0
                for j in idx[m][1]:idx[m][end]   
                    set1 += summed_strengths[a, j] * (system_geometry.beta[i, j] * system_geometry.sine_angle_panels[i, j] - log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.cos_angle_panels[i, j])
                    set2 += system_geometry.beta[i, j] * system_geometry.cos_angle_panels[i, j] + log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.sine_angle_panels[i, j]
                end
                # Swap indices to make AoA (a) the row and panel index (i) the column
                tangential_velocities[m][a, i - idx[m][1] + 1] = method.V_inf * (panel_geometry[m].cosine_vector[i] * cosd(flow_angles[a]) + panel_geometry[m].sine_vector[i] * sind(flow_angles[a])) + (set1 / (2 * π)) + (summed_strengths[a, end] / (2 * π)) * set2
                pressure_coefficient[m][a, i - idx[m][1] + 1] = 1.0 - (tangential_velocities[m][a, i - idx[m][1] + 1] / method.V_inf)^2
            end
        end
    end

    # @show tangential_velocities

    return (; tangential_velocities, pressure_coefficient)
end

# Loop bodies, angles of attack, and panels