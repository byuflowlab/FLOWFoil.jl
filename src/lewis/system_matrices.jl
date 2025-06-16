"""
    generate_system_matrices(method::Lewis, panel_geometry, system_geometry)

Generate the system matrices (coefficient matrix `A` and right-hand side vector `b`)
for the Lewis method given panel geometry and system geometry.

# Arguments
- `method::Lewis`: Marker type indicating the Lewis method.
- `panel_geometry`: A single panel geometry object.
- `system_geometry`: System geometry NamedTuple containing relative panel positions.

# Returns
- NamedTuple with fields:
  - `A::Matrix{Float}`: Coefficient matrix of the linear system.
  - `b::Vector{Float}`: Right-hand side vector for boundary conditions.
"""
function generate_system_matrices(method::Lewis, panel_geometry, system_geometry)
    return generate_system_matrices(method, [panel_geometry], system_geometry)
end

"""
    generate_system_matrices(method::Lewis, panel_geometry::AbstractVector, system_geometry)

Generate the system matrices (coefficient matrix `A` and right-hand side vector `b`)
for multiple bodies using the Lewis method.

# Arguments
- `method::Lewis`: Marker type indicating the Lewis method.
- `panel_geometry::Vector{Panel}`: Vector of panel geometry objects for each body.
- `system_geometry`: System geometry NamedTuple containing relative panel positions.

# Returns
- NamedTuple with fields:
  - `A::Matrix{Float}`: Coefficient matrix of the linear system.
  - `b::Vector{Float}`: Right-hand side vector for boundary conditions.
"""
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
    assemble_ring_vortex_matrix(body_of_revolution, panel_geometry, system_geometry)

Assembles the coefficient matrix `A` of influence coefficients for the
Lewis ring vortex method, applying back diagonal corrections and Kutta conditions.

# Arguments
- `body_of_revolution::Vector{Bool}`: Flags indicating if bodies are bodies of revolution.
- `panel_geometry::Vector{Panel}`: Vector of panel geometry objects.
- `system_geometry`: NamedTuple containing system geometry.

# Returns
- `A::Matrix{Float}`: The assembled influence coefficient matrix.
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
    get_kutta_indices(body_of_revolution, system_geometry)

Identify the panel indices for each body where Kutta conditions apply.
For bodies that are NOT bodies of revolution, returns leading and trailing edge panel indices.

# Arguments
- `body_of_revolution::Vector{Bool}`: Flags indicating if bodies are bodies of revolution.
- `system_geometry`: NamedTuple containing system geometry.

# Returns
- `kutta_idxs::Matrix{Int}`: Nx2 matrix where each row contains indices of
  leading and trailing edge panels for bodies requiring Kutta condition.
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
    apply_back_diagonal_correction!(amat, panel_geometry, idx, m2p)

Apply back diagonal correction to the coefficient matrix `amat` for bodies
that are not bodies of revolution. This modifies the matrix in-place.

# Arguments
- `amat::Matrix{Float}`: Coefficient matrix to modify.
- `panel_geometry::Panel`: Panel geometry object for the current body.
- `idx::UnitRange`: Index range for panels of the current body.
- `m2p::Vector{Int}`: Mapping from mesh indices to panel indices.

# Returns
- `nothing`: modifies `amat` in place.
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

"""
    assemble_ring_vortex_matrix_raw(body_of_revolution, panel_geometry, system_geometry)

Compute the raw influence coefficient matrix `A` before corrections for all bodies and panels.

# Arguments
- `body_of_revolution::Vector{Bool}`: Flags indicating if bodies are bodies of revolution.
- `panel_geometry::Vector{Panel}`: Vector of panel geometry objects.
- `system_geometry`: NamedTuple containing system geometry.

# Returns
- `amat::Matrix{Float}`: Raw influence coefficient matrix.
"""
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

Assemble the boundary condition vector `b` for the Lewis method, applying Kutta condition subtractions.

# Arguments
- `body_of_revolution::Vector{Bool}`: Flags indicating if bodies are bodies of revolution.
- `panel_geometry::Vector{Panel}`: Vector of panel geometry objects.
- `system_geometry`: NamedTuple containing system geometry.

# Returns
- `b::Vector{Float}`: Boundary condition vector with Kutta conditions applied.
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

"""
    assemble_ring_boundary_conditions_raw(body_of_revolution, panel_geometry, system_geometry)

Compute the raw boundary condition vector before applying Kutta condition subtractions.

# Arguments
- `body_of_revolution::Vector{Bool}`: Flags indicating if bodies are bodies of revolution.
- `panel_geometry::Vector{Panel}`: Vector of panel geometry objects.
- `system_geometry`: NamedTuple containing system geometry.

# Returns
- `bc::Vector{Float}`: Raw boundary condition vector.
"""
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
