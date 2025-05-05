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
    panelidx = system_geometry.mesh2panel

    # chord length
    chord = system_geometry.chord_length

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

    vs = [
        zeros(idx[m][end] - idx[m][1] + 1, naoa) for m in 1:nbodies
    ]
    cp = [
        zeros(idx[m][end] - idx[m][1] + 1, naoa) for m in 1:nbodies
    ]
    cl = zeros(naoa, nbodies)
    cd = zeros(naoa, nbodies)
    cm = zeros(naoa, nbodies)
    
    summed_strengths = zeros((maximum(idx[m][end] for m in 1:nbodies) + 1), naoa)

    for m in 1:nbodies  # Loop through number of bodies
        panel_vector = panel_geometry[m].panel_vectors
        for a in 1:naoa  # Loop through the different angles of attack         
            strengths_copy = copy(strengths)
            strengths_copy[:, 1] *= cosd(flow_angles[a])
            strengths_copy[:, 2] *= sind(flow_angles[a])
            summed_strengths[:, a] = sum(strengths_copy; dims=2) * 2 * pi * method.V_inf
            for i in idx[m][1]:idx[m][end]  # Loop through the panels
                set1 = 0.0
                set2 = 0.0
                for j in idx[m][1]:idx[m][end]   
                    set1 += summed_strengths[j, a] * (system_geometry.beta[i, j] * system_geometry.sine_angle_panels[i, j] - log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.cos_angle_panels[i, j])
                    set2 += system_geometry.beta[i, j] * system_geometry.cos_angle_panels[i, j] + log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.sine_angle_panels[i, j]
                end
                # Swap indices to make AoA (a) the row and panel index (i) the column
                vs[m][i - idx[m][1] + 1, a] = method.V_inf * (panel_geometry[m].cosine_vector[i] * cosd(flow_angles[a]) + panel_geometry[m].sine_vector[i] * sind(flow_angles[a])) + (set1 / (2 * π)) + (summed_strengths[end, a] / (2 * π)) * set2
                cp[m][i - idx[m][1] + 1, a] = 1.0 - (vs[m][i - idx[m][1] + 1, a])^2
            end

            ### --- Calculate Lift Coefficient --- ###
            cl[a, m] =
                sum([
                    cp[m][i, a] * (
                        -sind(flow_angles[a]) * panel_vector[i, 2] -
                        cosd(flow_angles[a]) * panel_vector[i, 1]
                    ) for i in panelidx[idx[m]]
                ]) / chord
        end
    end
    if nbodies == 1
        #if it is a single body, this reduces the need to use the body index
        vs_new = zeros(idx[1][end]-idx[1][1]+1, naoa)
        cp_new = zeros(idx[1][end]-idx[1][1]+1, naoa)
        cl_new = Array{TF}(undef, naoa) .= 0.0
        cd_new = similar(cl_new) .= 0.0
        cm_new = similar(cl_new) .= 0.0

        vs_new[:,:] = vs[1][:,:]
        cp_new[:,:] = cp[1][:,:]
        for i = 1:naoa
            cl_new[i] = cl[i,1]
            cd_new[i] = cd[i,1]
            cm_new[i] = cm[i,1]
        end

        vs = vs_new
        cp = cp_new
        cl = cl_new
        cd = cd_new
        cm = cm_new
    end
    return (; vs, cp, cl, cd, cm)
end

# Loop bodies, angles of attack, and panels