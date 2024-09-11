function generate_system_matrices(
    p::AxisymmetricProblem, panels::TP, mesh
) where {TP<:Panel}
    return generate_system_matrices(p, [panels], mesh)
end

function generate_system_matrices(p::AxisymmetricProblem, panels, mesh)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_ring_influence_matrix(p.singularity, p.body_of_revolution, panels, mesh)

    # Get boundary conditions (b, right hand side)
    b = assemble_ring_boundary_conditions(p.boundary, p.body_of_revolution, panels, mesh)

    return InviscidSystem(A, b, mesh.panel_indices)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#

"""
    assemble_ring_influence_matrix(v::Vortex, body_of_revolution, panels, mesh)

Assembles the "A" matrix (left hand side coefficient matrix).

**Arguments:**
- `s::Singularity` : The singularity type used.
- `body_of_revolution::Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

**Returns:**
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_ring_influence_matrix(::Singularity, body_of_revolution, panels, mesh) end

function assemble_ring_influence_matrix(v::Vortex, body_of_revolution, panels, mesh)
    return assemble_ring_vortex_matrix(v.order, body_of_revolution, panels, mesh)
end

"""
    assemble_ring_vortex_matrix(::Order, body_of_revolution, panels, mesh)

Assembles the coefficient matrix for a given order of singularity.

**Arguments:**
- `o::Order` : The order of singularity used.
- `body_of_revolution::Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

**Returns:**
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_ring_vortex_matrix(::Order, body_of_revolution, panels, mesh) end

function get_kutta_indices(body_of_revolution, mesh)

    # Count number of bodies requiring a Kutta Condition
    nk = count(br -> br == false, body_of_revolution)
    kutta_count = 1
    kutta_idxs = zeros(Int, nk, 2)

    for m in findall(m -> m == false, body_of_revolution)
        ### --- GetKutta Condition Indices --- ###
        kutta_idxs[kutta_count, :] = [mesh.panel_indices[m][1]; mesh.panel_indices[m][end]]
        kutta_count += 1
    end

    return kutta_idxs
end

function assemble_ring_vortex_matrix(::Constant, body_of_revolution, panels, mesh)
    amat = assemble_ring_vortex_matrix_raw(Constant(), body_of_revolution, panels, mesh)

    for m in findall(m -> m == false, body_of_revolution)
        apply_back_diagonal_correction!(
            amat, panels[m], mesh.panel_indices[m], mesh.mesh2panel
        )
    end

    kutta_idxs = get_kutta_indices(body_of_revolution, mesh)

    ## -- Apply Kutta Condition Subtractions -- ##
    for i in 1:length(kutta_idxs[:, 1])
        amat[kutta_idxs[i, 1], :] .-= amat[kutta_idxs[i, 2], :]
        amat[:, kutta_idxs[i, 1]] .-= amat[:, kutta_idxs[i, 2]]
    end

    return amat[1:end .∉ [kutta_idxs[:, 2]], 1:end .∉ [kutta_idxs[:, 2]]]
end

function apply_back_diagonal_correction!(amat, panels, idx, m2p)
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
                # sum += amat[j, i] * panels[m].panel_length[mesh2panel[j]]
                sum += amat[j, i] * panels.panel_length[m2p[j]]
            end
        end
        # dmagj = panels[m].panel_length[m2p[jidx]]
        dmagj = panels.panel_length[m2p[jidx]]
        amat[jidx, i] = -sum / dmagj
    end
    # end

    return nothing
end

function assemble_ring_vortex_matrix_raw(::Constant, body_of_revolution, panels, mesh)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies
    mesh2panel = mesh.mesh2panel

    # initialize coefficient matrix
    TF = eltype(mesh.m)
    # amat = zeros(TF, (N + nk, N + nk))
    amat = zeros(TF, (N, N))

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panels --- ###
            for i in idx[m]
                for j in idx[n]

                    ### --- Calculate influence coefficient --- ###
                    amat[i, j] = calculate_ring_vortex_influence(
                        Constant(), panels[m], panels[n], mesh, i, j
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
    assemble_ring_boundary_conditions(::BoundaryCondition, body_of_revolution, panels, mesh)

Assemble boundary condition vector.

**Arguments:**
- `bc::BoundaryCondition` : The type of boundary condition to be used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `body_of_revolution::Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

**Returns**
 - `b::Matrix{Float}` : Boundary condition matrix
"""
function assemble_ring_boundary_conditions(
    ::BoundaryCondition, body_of_revolution, panels, mesh
) end

function assemble_ring_boundary_conditions(::Dirichlet, body_of_revolution, panels, mesh)
    bc = assemble_ring_boundary_conditions_raw(
        Dirichlet(), body_of_revolution, panels, mesh
    )

    kutta_idxs = get_kutta_indices(body_of_revolution, mesh)

    ## -- Apply Kutta Condition Subtractions -- ##
    for i in 1:length(kutta_idxs[:, 1])
        bc[kutta_idxs[i, 1]] -= bc[kutta_idxs[i, 2]]
    end

    return bc[1:end .∉ [kutta_idxs[:, 2]]]
end

function assemble_ring_boundary_conditions_raw(
    ::Dirichlet, body_of_revolution, panels, mesh
)

    ### --- SETUP --- ###

    # # Count number of bodies requiring a Kutta Condition
    # nk = count(br -> br == false, body_of_revolution)
    # kutta_idxs = zeros(Int, nk, 2)
    # kutta_count = 1

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies
    mesh2panel = mesh.mesh2panel

    # initialize boundary condition array
    TF = eltype(mesh.m)
    # bc = zeros(TF, N + nk)
    bc = zeros(TF, N)

    ### --- Loop through bodies --- ###
    for m in 1:nbodies

        # generate portion of boundary condition array associated with mth body
        bc[idx[m], 1] = [-cos(panels[m].panel_angle[mesh2panel[i]]) for i in idx[m]]
        # if !body_of_revolution[m]
        #     kutta_idxs[kutta_count, :] = [idx[m][1]; idx[m][end]]
        #     kutta_count += 1
        # end
    end

    return bc

    # ## -- Apply Kutta Condition Subtractions -- ##
    # for i in 1:nk
    #     bc[kutta_idxs[i, 1]] -= bc[kutta_idxs[i, 2]]
    # end

    # return bc[1:end .∉ [kutta_idxs[:, 2]]]
end
