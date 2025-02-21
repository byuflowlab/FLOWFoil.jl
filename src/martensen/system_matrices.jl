function generate_system_matrices(method::Martensen, panel_geometry, system_geometry)
    return generate_system_matrices(method, [panel_geometry], system_geometry)
end

function generate_system_matrices(
    method::Martensen, panel_geometry::AbstractVector, system_geometry
)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_periodic_vortex_matrix(
        panel_geometry,
        system_geometry,
        (;
            method.pitch,
            method.stagger,
            method.cascade,
            method.transition_value,
            method.transition_hardness,
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

# Arguments:
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

# Returns:
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

                    ### --- Calculate influence coefficient --- ###

                    # - Self-induced coefficient - #
                    if i == j
                        amat[i, j] = calculate_periodic_self_vortex_influence(
                            panel_geometry[m], i, cascade_parameters.curvature_correction
                        )

                    else
                        # - Non self-induced coefficient - #

                        # Calculate planar influence coefficient
                        K_planar = calculate_planar_vortex_influence(
                            panel_geometry[m],
                            panel_geometry[n],
                            system_geometry,
                            i,
                            j,
                            cascade_parameters,
                        )

                        # if desired to include high solidity cascade effects:
                        if cascade_parameters.cascade

                            # Calculate periodic influence coefficient
                            K_periodic = calculate_periodic_vortex_influence(
                                panel_geometry[m],
                                panel_geometry[n],
                                system_geometry,
                                i,
                                j,
                                cascade_parameters,
                            )

                            # blend coefficients based on transition value
                            amat[i, j] = FLOWMath.sigmoid_blend(
                                K_periodic,
                                K_planar,
                                system_geometry.pitch_to_chord,
                                cascade_parameters.transition_value,
                                cascade_parameters.transition_hardness,
                            )
                        else

                            # if only planar analysis is desired, skip periodic stuff
                            amat[i, j] = K_planar
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

# Arguments:
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
