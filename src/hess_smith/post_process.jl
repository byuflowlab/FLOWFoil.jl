function post_process(method::HessSmith, panel_geometry, system_geometry, strengths; npanels=80)
    return post_process(
        method, [panel_geometry], system_geometry, strengths; npanels=npanels
    )
end

"""

Calculate tangential velocity and surface pressures

"""
function post_process(
    method::HessSmith, panel_geometry::AbstractVector, system_geometry, strengths; npanels=80
)

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    println(idx)
    nbodies = system_geometry.nbodies

    # - Initialize Outputs - #
    TF = eltype(system_geometry.k2)

    tangential_velocities = [zeros(idx[m][end]-idx[m][1]+1) for m in 1:nbodies]
    surface_pressures = [zeros(idx[m][end]-idx[m][1]+1) for m in 1:nbodies]

    for m in 1:nbodies
        # - Extract surface velocity - #
        #= Note that we here assume that we are using the subtractive method for the kutta conditions, requiring us to recover the strengths value for the last panel (trailing edge upper side).  We also assume here that the indexing starts at the lower side trailing edge and proceeds clockwise back to the upper side trailing edge.  Otherwise, not only will the solver not have worked, there will also be an indexing error here
        =#

        for j in eachindex(strengths[1:end-1])    
            set1 += strengths[j] * (system_geometry.beta[i, j] * system_geometry.sine_angle_panels[i, j] - log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.cos_angle_panels[i, j])
            set2 += system_geometry.beta[i, j] * system_geometry.cos_angle_panels[i, j] + log(system_geometry.r_influence[i, j+1] / system_geometry.r_influence[i, j]) * system_geometry.sine_angle_panels[i, j]
        end
        tangential_velocity[m][:] = V_inf * (panel_geometry.cosine_vector[i] * cos(panel_geometry.AoA) + panel_geometry.sine_vector[i] * sin(panel_geometry.AoA)) + (set1 / (2 * π)) + (strengths[end] / (2 * π)) * set2

        # - Calculate surface pressure - #
        surface_pressures[m][:] = 1.0 .- (tangential_velocities[m][:]) .^ 2
    end

    return (; tangential_velocities, surface_pressures)
end