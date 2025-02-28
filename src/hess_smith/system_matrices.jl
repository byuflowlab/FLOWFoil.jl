function generate_system_matrices(method::HessSmith, panel_geometry, system_geometry)
    return generate_system_matrices(method, [panel_geometry], system_geometry)
end

function generate_system_matrices(
    method::HessSmith, panel_geometry::AbstractVector, system_geometry
)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_vortex_matrix(
        panel_geometry,
        system_geometry
        )

    # Get boundary conditions (b, right hand side)
    b = assemble_b_matrix(
        panel_geometry,
        system_geometry
        )

    return (; A, b)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#
"""
    assemble_vortex_matrix(panel_geometry, system_geometry)

Assembles the coefficient matrix for a given order of singularity.

# Arguments:
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

# Returns:
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_vortex_matrix(
    panel_geometry, system_geometry
)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    N = idx[end][end]
    nbodies = system_geometry.nbodies

    # initialize coefficient matrix
    TF = eltype(system_geometry.r_x)
    amat = zeros(TF, (N, N))

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panel_geometry --- ###
            for i in idx[m]
                k1 = system_geometry[m].beta[1, i] * system_geometry[m].sine_angle_panels[1, i] - log(system_geometry[m].r_influence[1, i+1] / system_geometry[m].r_influence[1, i]) * system_geometry[m].cosine_angle_panels[1, i]
                kn = system_geometry[m].beta[end, i] * system_geometry[m].sine_angle_panels[end, i] - log(system_geometry[m].r_influence[end, i+1] / system_geometry[m].r_influence[end, i]) * system_geometry[m].cosine_angle_panels[end, i]
                A[end, i] = k1 + kn
                A[end, end] = system_geometry[m].beta[end, i] * system_geometry[m].cosine_angle_panels[end, i] + log(system_geometry[m].r_influence[end, i+1] / system_geometry[m].r_influence[end, i]) * system_geometry[m].sine_angle_panels[end, i] + system_geometry[m].beta[1, i] * system_geometry[m].cosine_angle_panels[1, i] + log(system_geometry[m].r_influence[1, i+1] / system_geometry[m].r_influence[1, i]) * system_geometry[m].sine_angle_panels[1, i]
                
                for j in idx[n]
                    A[i, j] = log(system_geometry[m].r_influence[i, j+1] / system_geometry[m].r_influence[i, j]) * system_geometry[m].sine_angle_panels[i, j] + system_geometry[m].beta[i, j] * system_geometry[m].cosine_angle_panels[i, j]
                    A[i, end] += log(system_geometry[m].r_influence[i, j+1] / system_geometry[m].r_influence[i, j]) * system_geometry[m].cosine_angle_panels[i, j] - system_geometry[m].beta[i, j] * system_geometry[m].sine_angle_panels[i, j] 
                end
            end
        end
    end

    return amat
end

"""
"""
function assemble_b_matrix(panel_geometry, system_geometry)
    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    N = idx[end][end]
    nbodies = system_geometry.nbodies

    # initialize coefficient matrix
    TF = eltype(system_geometry.r_x)
    bmat = zeros(TF, N)
    V_inf = 1.0

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panel_geometry --- ###
            for i in idx[m]
                bmat[i] = 2 * π * V_inf * (panel_geometry[m].sine_vector[i] * cos(panel_geometry.AoA) - panel_geometry[m].cosine_vector[i] * sin(panel_geometry.AoA)) 
            end
            bmat[end] = -2 * π * V_inf * ((panel_geometry[m].cosine_vector[1] * cos(panel_geometry.AoA) + panel_geometry[m].sine_vectorl[1] * sin(panel_geometry.AoA)) + (panel_geometry[m].cosine_vector[end] * cos(panel_geometry.AoA) + panel_geometry[m].sine_vector[end] * sin(panel_geometry.AoA)))
        end
    end

    return bmat
end
### I still need to figure out where we input V_inf and panel_geometry.AoA, though to be honest, V_inf doesn't matter since it gets non-dimensionalized later.
### AoA is currently input as an argument into the generate_panel_geometry function