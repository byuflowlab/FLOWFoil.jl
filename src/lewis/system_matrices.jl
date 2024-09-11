function generate_system_matrices(method::Lewis, panel_geometry, system_geometry)
    return generate_system_matrices(method, [panel_geometry], system_geometry)
end

function generate_system_matrices(
    method::Lewis, panel_geometry::AbstractVector, system_geometry
)

    # Get coeffiecient matrix (A, left hand side)
    A = assemble_ring_vortex_matrix(
        method.body_of_revolution, panel_geometry, system_geometry
    )

    # Get boundary conditions (b, right hand side)
    b = assemble_right_hand_side(method.body_of_revolution, panel_geometry, system_geometry)

    return (; A, b)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#

"""
    assemble_ring_influence_matrix(body_of_revolution, panel_geometry, system_geometry)

Assembles the "A" matrix (left hand side coefficient matrix).

**Arguments:**
- `body_of_revolution::Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

**Returns:**
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_ring_vortex_matrix(body_of_revolution, panel_geometry, system_geometry)
    amat = assemble_ring_vortex_matrix_raw(
        body_of_revolution, panel_geometry, system_geometry
    )

    for m in findall(m -> m == false, body_of_revolution)
        apply_back_diagonal_correction!(
            amat,
            panel_geometry[m],
            system_geometry.panel_indices[m],
            system_geometry.mesh2panel,
        )
    end

    kutta_idxs = get_kutta_indices(body_of_revolution, system_geometry)

    ## -- Apply Kutta Condition Subtractions -- ##
    for i in 1:length(kutta_idxs[:, 1])
        amat[kutta_idxs[i, 1], :] .-= amat[kutta_idxs[i, 2], :]
        amat[:, kutta_idxs[i, 1]] .-= amat[:, kutta_idxs[i, 2]]
    end

    return amat[1:end .∉ [kutta_idxs[:, 2]], 1:end .∉ [kutta_idxs[:, 2]]]
end

"""
"""
function get_kutta_indices(body_of_revolution, system_geometry)

    # Count number of bodies requiring a Kutta Condition
    nk = count(br -> br == false, body_of_revolution)
    kutta_count = 1
    kutta_idxs = zeros(Int, nk, 2)

    for m in findall(m -> m == false, body_of_revolution)
        ### --- GetKutta Condition Indices --- ###
        kutta_idxs[kutta_count, :] = [
            system_geometry.panel_indices[m][1]
            system_geometry.panel_indices[m][end]
        ]
        kutta_count += 1
    end

    return kutta_idxs
end

"""
"""
function apply_back_diagonal_correction!(amat, panel_geometry, idx, m2p)
    # if m == n && !body_of_revolution[m]

    ### --- Apply Back Diagonal Correction --- ###
    # for i in idx[m]
    for i in idx
        sum = 0.0
        # jidx = idx[m][end] + 1 - i
        jidx = idx[end] + 1 - i
        # for j in idx[m]
        for j in idx
            if j != jidx
                # sum += amat[j, i] * panel_geometry[m].panel_length[mesh2panel[j]]
                sum += amat[j, i] * panel_geometry.panel_length[m2p[j]]
            end
        end
        # dmagj = panel_geometry[m].panel_length[m2p[jidx]]
        dmagj = panel_geometry.panel_length[m2p[jidx]]
        amat[jidx, i] = -sum / dmagj
    end
    # end

    return nothing
end

function assemble_ring_vortex_matrix_raw(
    body_of_revolution, panel_geometry, system_geometry
)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    N = idx[end][end]
    nbodies = system_geometry.nbodies
    mesh2panel = system_geometry.mesh2panel

    # initialize coefficient matrix
    TF = eltype(system_geometry.k2)
    # amat = zeros(TF, (N + nk, N + nk))
    amat = zeros(TF, (N, N))

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panel_geometry --- ###
            for i in idx[m]
                for j in idx[n]

                    ### --- Calculate influence coefficient --- ###
                    amat[i, j] = calculate_ring_vortex_influence(
                        panel_geometry[m], panel_geometry[n], system_geometry, i, j
                    )
                end
            end
        end
    end

    return amat

    # ## -- Apply Kutta Condition Subtractions -- ##
    # for i in 1:nk
    #     amat[kutta_idxs[i, 1], :] .-= amat[kutta_idxs[i, 2], :]
    #     amat[:, kutta_idxs[i, 1]] .-= amat[:, kutta_idxs[i, 2]]
    # end

    # return amat[1:end .∉ [kutta_idxs[:, 2]], 1:end .∉ [kutta_idxs[:, 2]]]
end

#---------------------------------#
#    BOUNDARY CONDITION MATRIX    #
#---------------------------------#

"""
    assemble_right_hand_side(body_of_revolution, panel_geometry, system_geometry)

Assemble boundary condition vector.

**Arguments:**
- `panel_geometry::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `body_of_revolution::Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `system_geometry::system_geometry` : The system_geometry object containing relative geometry for the influence coefficient calculations.

**Returns**
 - `b::Matrix{Float}` : Boundary condition matrix
"""
function assemble_right_hand_side(body_of_revolution, panel_geometry, system_geometry)
    bc = assemble_ring_boundary_conditions_raw(
        body_of_revolution, panel_geometry, system_geometry
    )

    kutta_idxs = get_kutta_indices(body_of_revolution, system_geometry)

    ## -- Apply Kutta Condition Subtractions -- ##
    for i in 1:length(kutta_idxs[:, 1])
        bc[kutta_idxs[i, 1]] -= bc[kutta_idxs[i, 2]]
    end

    return bc[1:end .∉ [kutta_idxs[:, 2]]]
end

function assemble_ring_boundary_conditions_raw(
    body_of_revolution, panel_geometry, system_geometry
)

    ### --- SETUP --- ###

    # # Count number of bodies requiring a Kutta Condition
    # nk = count(br -> br == false, body_of_revolution)
    # kutta_idxs = zeros(Int, nk, 2)
    # kutta_count = 1

    # - Rename for Convenience - #
    idx = system_geometry.panel_indices
    N = idx[end][end]
    nbodies = system_geometry.nbodies
    mesh2panel = system_geometry.mesh2panel

    # initialize boundary condition array
    TF = eltype(system_geometry.k2)
    # bc = zeros(TF, N + nk)
    bc = zeros(TF, N)

    ### --- Loop through bodies --- ###
    for m in 1:nbodies

        # generate portion of boundary condition array associated with mth body
        bc[idx[m], 1] = [-cos(panel_geometry[m].panel_angle[mesh2panel[i]]) for i in idx[m]]
        # if !body_of_revolution[m]
        #     kutta_idxs[kutta_count, :] = [idx[m][1]; idx[m][end]]
        #     kutta_count += 1
        # end
    end

    return bc
end
