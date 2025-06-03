"""
    post_process(method::HessSmith, panel_geometry, system_geometry, strengths, flow_angles)

Wrapper function for computing post-processing aerodynamic results when only a single body
is used in the panel method formulation. Converts a single-body input into a vector format expected by the full
multi-body `post_process` method, ensuring uniform handling for downstream calculations.

# Arguments
- `method::HessSmith`: The solver type used to dispatch Hess-Smith panel method logic.
- `panel_geometry`: The panel geometry object for a single body (not a vector).
- `system_geometry`: Geometry-related precomputed fields needed for influence calculations.
- `strengths::Matrix`: Matrix of vortex strengths per panel.
- `flow_angles::Vector{<:Real}`: Array of freestream angles of attack (in degrees).

# Returns
- An `InviscidOutputs` object containing:
  - `vs`: Tangential velocities.
  - `cp`: Pressure coefficients.
  - `cl`: Lift coefficients.
  - `cd`: Drag coefficients (zeroed for inviscid).
  - `cm`: Moment coefficients (zeroed for inviscid).
"""
function post_process(
    method::HessSmith, panel_geometry, system_geometry, strengths, flow_angles
)
    return post_process(method, [panel_geometry], system_geometry, strengths, flow_angles)
end

"""
    post_process(::HessSmith, panel_geometry, system_geometry, strengths, flow_angles)

Computes tangential velocities, pressure coefficients, and aerodynamic coefficients
(lift, drag, and moment) for a given configuration of panels and freestream conditions
using the Hess-Smith panel method.

# Arguments
- `::HessSmith`: Marker type for dispatching the Hess-Smith solver method.
- `panel_geometry::AbstractVector`: Vector of panel geometry objects, each representing an airfoil or body.
- `system_geometry`: Precomputed system geometry structure (distances, angles, etc.) used for influence calculations.
- `strengths::Matrix`: Matrix of vortex strengths for each panel and each vector component.
- `flow_angles::Vector{<:Real}`: Freestream angles of attack (in degrees) to evaluate aerodynamic performance.

# Returns
- If `nbodies == 1`: An `InviscidOutputs` struct containing:
  - `vs`: Matrix of tangential velocities for each panel and angle of attack.
  - `cp`: Matrix of surface pressure coefficients.
  - `cl`: Vector of lift coefficients.
  - `cd`: Vector of drag coefficients (currently zero for inviscid flow).
  - `cm`: Vector of moment coefficients (currently zero for inviscid flow).

- If `nbodies > 1`: An `InviscidOutputs` struct containing:
  - `vs::Vector{Matrix}`: List of tangential velocity matrices for each body.
  - `cp::Vector{Matrix}`: List of pressure coefficient matrices for each body.
  - `cl::Matrix`: Matrix of lift coefficients `[n_aoa, n_bodies]`.
  - `cd::Matrix`: Matrix of drag coefficients (zero).
  - `cm::Matrix`: Matrix of moment coefficients (zero).
"""
function post_process(
    method::HessSmith,
    panel_geometry::AbstractVector,
    system_geometry,
    strengths,
    flow_angles,
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
    # #     summed_strengths[:, a] = sum(strengths; dims=2) * 2.0* pi * method.V_inf
    # # end
    # for a in 1:naoa
    #     strengths_copy = copy(strengths)
    #     strengths_copy[:, 1] *= cosd(flow_angles[a])
    #     strengths_copy[:, 2] *= sind(flow_angles[a])
    #     summed_strengths[:, a] = sum(strengths_copy; dims=2) * 2.0* pi * method.V_inf
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
    #             tangential_velocities[m][i - idx[m][1] + 1, a] = method.V_inf * (panel_geometry[m].cosine_vector[i] * cosd(flow_angles[a]) + panel_geometry[m].sine_vector[i] * sind(flow_angles[a])) + (set1 / (2.0* pi)) + (summed_strengths[end, a] / (2.0* pi)) * set2
    #             # surface_pressures[m][i - idx[m][1] + 1, a] = 1.0 - tangential_velocities[m][i - idx[m][1] + 1, a] ^ 2
    #             pressure_coefficient[m][i-idx[m][1] + 1, a] = 1.0 - (tangential_velocities[m][i - idx[m][1] + 1, a] / method.V_inf)^2
    #         end
    #     end
    # end

    vs = [zeros(idx[m][end] - idx[m][1] + 1, naoa) for m in 1:nbodies]
    cp = [zeros(idx[m][end] - idx[m][1] + 1, naoa) for m in 1:nbodies]
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
            summed_strengths[:, a] = sum(strengths_copy; dims=2) * 2.0 * pi * method.V_inf
            for i in idx[m][1]:idx[m][end]  # Loop through the panels
                set1 = 0.0
                set2 = 0.0
                for j in idx[m][1]:idx[m][end]
                    set1 +=
                        summed_strengths[j, a] * (
                            system_geometry.beta[i, j] *
                            system_geometry.sine_angle_panels[i, j] -
                            log(
                                system_geometry.r_influence[i, j + 1] /
                                system_geometry.r_influence[i, j],
                            ) * system_geometry.cos_angle_panels[i, j]
                        )
                    set2 +=
                        system_geometry.beta[i, j] *
                        system_geometry.cos_angle_panels[i, j] +
                        log(
                            system_geometry.r_influence[i, j + 1] /
                            system_geometry.r_influence[i, j],
                        ) * system_geometry.sine_angle_panels[i, j]
                end
                # Swap indices to make AoA (a) the row and panel index (i) the column
                vs[m][i - idx[m][1] + 1, a] =
                    method.V_inf * (
                        panel_geometry[m].cosine_vector[i] * cosd(flow_angles[a]) +
                        panel_geometry[m].sine_vector[i] * sind(flow_angles[a])
                    ) +
                    (set1 / (2.0 * pi)) +
                    (summed_strengths[end, a] / (2.0 * pi)) * set2
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
        return InviscidOutputs(vs[1], cp[1], cl, cd, cm)
    else
        return InviscidOutputs(vs, cp, cl, cd, cm)
    end
end

# Loop bodies, angles of attack, and panels
