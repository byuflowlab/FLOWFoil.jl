"""
    generate_system_matrices(method::HessSmith, panel_geometry, system_geometry)

Wrapper function to generate the Hess-Smith system matrices for a single body. Internally wraps
the panel geometry into a vector and delegates to the general multi-body function.

# Arguments
- `method::HessSmith`: The solver object or configuration struct for the Hess-Smith panel method.
- `panel_geometry`: A named tuple representing a single body's panel geometry (e.g., from `generate_panel_geometry`).
- `system_geometry`: A named tuple containing geometric data between panels (e.g., from `generate_system_geometry`).

# Returns
- A named tuple `(; A, b)` where:
  - `A::Matrix`: The system matrix representing the vortex influence coefficients (LHS of linear system).
  - `b::Vector`: The boundary condition vector (RHS of linear system).
"""
function generate_system_matrices(method::HessSmith, panel_geometry, system_geometry)
    return generate_system_matrices(method, [panel_geometry], system_geometry)
end

"""
    generate_system_matrices(method::HessSmith, panel_geometry::AbstractVector, system_geometry)

Generates the system matrices for solving the Hess-Smith panel method. This includes the vortex
influence coefficient matrix (`A`) and the boundary condition vector (`b`), for one or more bodies.

# Arguments
- `method::HessSmith`: The solver object or configuration struct for the Hess-Smith panel method.
- `panel_geometry::AbstractVector`: A vector of named tuples representing the panel geometry of each body.
- `system_geometry`: A named tuple containing system-level geometric data (from `generate_system_geometry`).

# Returns
- A named tuple `(; A, b)` where:
  - `A::Matrix`: The system matrix of vortex influence coefficients.
  - `b::Vector`: The vector of boundary conditions imposed at the control points.
"""
function generate_system_matrices(
    method::HessSmith, panel_geometry::AbstractVector, system_geometry
)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_vortex_matrix(panel_geometry, system_geometry)

    # Get boundary conditions (b, right hand side)
    b = assemble_b_matrix(panel_geometry, system_geometry)

    return (; A, b)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#
"""
    assemble_vortex_matrix(panel_geometry, system_geometry)

Assembles the coefficient matrix.

# Arguments
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

# Returns
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_vortex_matrix(panel_geometry, system_geometry)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    N = idx[end][end]
    nbodies = system_geometry.nbodies

    # initialize coefficient matrix
    TF = eltype(system_geometry.r_x)
    amat = zeros(TF, (N + 1, N + 1))

    # Loop through system

    ### --- Loop through panel_geometry --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            for i in 1:N
                k1 =
                    system_geometry.beta[1, i] * system_geometry.sine_angle_panels[1, i] -
                    log(
                        system_geometry.r_influence[1, i + 1] /
                        system_geometry.r_influence[1, i],
                    ) * system_geometry.cos_angle_panels[1, i]
                kn =
                    system_geometry.beta[end, i] *
                    system_geometry.sine_angle_panels[end, i] -
                    log(
                        system_geometry.r_influence[end, i + 1] /
                        system_geometry.r_influence[end, i],
                    ) * system_geometry.cos_angle_panels[end, i]
                amat[end, i] = k1 + kn
                amat[end, end] =
                    system_geometry.beta[end, i] *
                    system_geometry.cos_angle_panels[end, i] +
                    log(
                        system_geometry.r_influence[end, i + 1] /
                        system_geometry.r_influence[end, i],
                    ) * system_geometry.sine_angle_panels[end, i] +
                    system_geometry.beta[1, i] * system_geometry.cos_angle_panels[1, i] +
                    log(
                        system_geometry.r_influence[1, i + 1] /
                        system_geometry.r_influence[1, i],
                    ) * system_geometry.sine_angle_panels[1, i]

                for j in 1:N
                    amat[i, j] =
                        log(
                            system_geometry.r_influence[i, j + 1] /
                            system_geometry.r_influence[i, j],
                        ) * system_geometry.sine_angle_panels[i, j] +
                        system_geometry.beta[i, j] * system_geometry.cos_angle_panels[i, j]
                    amat[i, end] +=
                        log(
                            system_geometry.r_influence[i, j + 1] /
                            system_geometry.r_influence[i, j],
                        ) * system_geometry.cos_angle_panels[i, j] -
                        system_geometry.beta[i, j] * system_geometry.sine_angle_panels[i, j]
                end
            end
        end
    end
    return amat
end

"""
    assemble_b_matrix(panel_geometry, system_geometry)

Assembles the boundary condition vector.

# Arguments
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

# Returns
- b::Vector{Float}` : The boundary condition vector for the linear system
"""
function assemble_b_matrix(panel_geometry, system_geometry)   ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    N = idx[end][end]
    nbodies = system_geometry.nbodies

    # initialize coefficient matrix
    TF = eltype(system_geometry.r_x)
    # bmat = [zeros(2) for _ in 1:(N+1)]
    bmat = zeros(TF, (N + 1, 2))

    # ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            N = idx[m][end]
            for i in 1:N
                panel_index_i = i
                if m > 1 #this is added to call the correct panel for the periodic vortex influence - which doesn't account for the multiple bodies
                    panel_index_i = i - panel_indices[m - 1][end]
                end
                bmat[i, :] = [
                    panel_geometry[m].sine_vector[panel_index_i]
                    -panel_geometry[m].cosine_vector[panel_index_i]
                ]
            end

            bmat[end, :] = [
                -(
                    panel_geometry[m].cosine_vector[1] +
                    panel_geometry[m].cosine_vector[end]
                )
                -(panel_geometry[m].sine_vector[1] + panel_geometry[m].sine_vector[end])
            ]
        end
    end

    return bmat
end
