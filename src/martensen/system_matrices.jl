"""
    generate_system_matrices(method::Martensen, panel_geometry, system_geometry)

Convenience method to generate system matrices for a single panel geometry body.

# Arguments
- `method::Martensen`: The Martensen method configuration object.
- `panel_geometry`: A single panel geometry object representing one body.
- `system_geometry`: Precomputed system geometry object for the panels.

# Returns
- Named tuple containing:
  - `A`: Coefficient matrix of the system.
  - `b`: Boundary condition vector (right-hand side).
"""
function generate_system_matrices(method::Martensen, panel_geometry, system_geometry)
    return generate_system_matrices(method, [panel_geometry], system_geometry)
end

"""
    generate_system_matrices(method::Martensen, panel_geometry::AbstractVector, system_geometry)

Generates the linear system matrices for multiple panel geometry bodies based on the Martensen method.

# Arguments
- `method::Martensen`: The Martensen method configuration object containing method parameters.
- `panel_geometry::AbstractVector`: Vector of panel geometry objects, one for each body.
- `system_geometry`: Precomputed system geometry object describing relative panel positions and distances.

# Returns
- Named tuple containing:
  - `A`: Coefficient matrix for the linear system.
  - `b`: Boundary condition vector (right-hand side).
"""
function generate_system_matrices(
    method::Martensen, panel_geometry::AbstractVector, system_geometry
)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_periodic_vortex_matrix(
        panel_geometry,
        system_geometry,
        (;
            method.solidity,
            method.stagger,
            method.cascade,
            method.transition_value,
            method.curvature_correction,
        ),
    )
    # Get boundary conditions (b, right hand side)
    b = assemble_periodic_right_hand_side(panel_geometry, system_geometry)
    return (; A, b)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#
"""
    assemble_periodic_vortex_matrix(panel_geometry, system_geometry, cascade_parameters)

Constructs the vortex influence coefficient matrix (A) for a system of panels, applying cascade effects
and Kutta conditions for each body.

# Arguments
- `panel_geometry`: Vector of panel geometry objects.
- `system_geometry`: System geometry containing panel relative distances and indexing.
- `cascade_parameters`: Named tuple containing parameters like `solidity`, `stagger`, `cascade`, 
  `transition_value`, and `curvature_correction`.

# Returns
- Matrix `amat` of influence coefficients with Kutta condition rows and columns appropriately modified and reduced.
"""
function assemble_periodic_vortex_matrix(
    panel_geometry, system_geometry, cascade_parameters
)
    amat = assemble_periodic_vortex_matrix_raw(
        panel_geometry, system_geometry, cascade_parameters
    )

    for m in 1:(system_geometry.nbodies)
        #NOTE: this back diagonal correction is defined in src/lewis/system_matrices.jl
        apply_back_diagonal_correction!(
            amat,
            panel_geometry[m],
            system_geometry.panel_indices[m],
            system_geometry.mesh2panel,
        )
    end

    kutta_idxs = get_kutta_indices(system_geometry)

    ## -- Apply Kutta Condition Subtractions -- ##
    for kid in eachrow(kutta_idxs)
        amat[kid[1], :] .-= amat[kid[2], :]
        amat[:, kid[1]] .-= amat[:, kid[2]]
    end
    return amat[1:end .∉ [kutta_idxs[:, 2]], 1:end .∉ [kutta_idxs[:, 2]]]
end

"""
    get_kutta_indices(system_geometry)

Computes the panel indices associated with the Kutta condition (leading and trailing edges) for each body.

# Arguments
- `system_geometry`: System geometry object with body and panel indexing information.

# Returns
- `kutta_idxs::Matrix{Int}`: A matrix of size (nbodies, 2), where each row contains the indices of the
  first and last panel for that body, representing Kutta condition panels.
"""
function get_kutta_indices(system_geometry)

    # Count number of bodies requiring a Kutta Condition
    nk = system_geometry.nbodies
    kutta_idxs = zeros(Int, nk, 2)

    for m in 1:nk
        ### --- GetKutta Condition Indices --- ###
        kutta_idxs[m, :] = [
            system_geometry.panel_indices[m][1]
            system_geometry.panel_indices[m][end]
        ]
    end

    return kutta_idxs
end

"""
    assemble_periodic_vortex_matrix_raw(panel_geometry, system_geometry)

Assembles the coefficient matrix for a given order of singularity.

# Arguments
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

# Returns
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_periodic_vortex_matrix_raw(
    panel_geometry, system_geometry, cascade_parameters
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
                for j in idx[n]
                    panel_index_i = i
                    panel_index_j = j
                    if m > 1 #this is added to call the correct panel for the periodic vortex influence - which doesn't account for the multiple bodies
                        panel_index_i = i - idx[m - 1][end]
                    end
                    if n > 1
                        panel_index_j = j - idx[n - 1][end]
                    end

                    ### --- Calculate influence coefficient --- ###

                    # - Self-induced coefficient - #
                    if panel_index_i == panel_index_j
                        # same for both cascade and non-cascade, just needed a unique function name
                        amat[i, j] = calculate_periodic_self_vortex_influence(
                            panel_geometry[m], panel_index_i, cascade_parameters.curvature_correction
                        )

                    else
                        # - Non self-induced coefficient - #
                        if cascade_parameters.cascade
                            # Calculate periodic influence coefficient
                            amat[i, j] = calculate_periodic_vortex_influence(
                                panel_geometry[m],
                                panel_geometry[n],
                                system_geometry,
                                panel_index_i,
                                panel_index_j,
                                cascade_parameters,
                            )

                        else
                            # Calculate planar influence coefficient
                            amat[i, j] = calculate_planar_vortex_influence(
                                panel_geometry[m],
                                panel_geometry[n],
                                system_geometry,
                                panel_index_i,
                                panel_index_j,
                                cascade_parameters,
                            )
                        end
                    end
                end
            end
        end
    end
    return amat
end

#---------------------------------#
#    BOUNDARY CONDITION MATRIX    #
#---------------------------------#

"""
    assemble_right_hand_side(panel_geometry, system_geometry)

Assemble boundary condition vector.

# Arguments
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

# Returns
 - `b::Matrix{Float}` : Boundary condition matrix
"""
function assemble_periodic_right_hand_side(panel_geometry, system_geometry)
    bc = assemble_periodic_boundary_conditions_raw(panel_geometry, system_geometry)

    kutta_idxs = get_kutta_indices(system_geometry)

    ## -- Apply Kutta Condition Subtractions -- ##
    for kid in eachrow(kutta_idxs)
        bc[kid[1], :] .-= bc[kid[2], :]
    end

    return bc[1:end .∉ [kutta_idxs[:, 2]], :]
end

"""
    assemble_periodic_boundary_conditions_raw(panel_geometry, system_geometry)

Generates the raw boundary condition matrix representing flow tangency on each panel without Kutta condition enforcement.

# Arguments
- `panel_geometry`: Vector of panel geometry objects.
- `system_geometry`: System geometry containing panel indexing and relative geometry.

# Returns
- `bc::Matrix{Float}`: Matrix where each row corresponds to a panel, containing the negative cosine and sine
  of the panel angle components representing flow normal velocity boundary conditions.
"""
function assemble_periodic_boundary_conditions_raw(panel_geometry, system_geometry)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    N = idx[end][end]
    nbodies = system_geometry.nbodies
    mesh2panel = system_geometry.mesh2panel

    # initialize boundary condition array
    TF = eltype(system_geometry.r_x)
    bc = zeros(TF, N, 2)

    ### --- Loop through bodies --- ###
    for m in 1:nbodies

        # generate portion of boundary condition array associated with mth body
        bc[idx[m], 1] = [-cos(panel_geometry[m].panel_angle[mesh2panel[i]]) for i in idx[m]]
        bc[idx[m], 2] = [-sin(panel_geometry[m].panel_angle[mesh2panel[i]]) for i in idx[m]]
    end

    return bc
end
