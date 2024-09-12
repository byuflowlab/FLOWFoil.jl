function generate_system_matrices(p::PeriodicProblem, panels::TP, mesh) where {TP<:Panel}
    return generate_system_matrices(p, [panels], mesh)
end

function generate_system_matrices(p::PeriodicProblem, panels, mesh)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_periodic_influence_matrix(p.singularity, panels, mesh)

    # Get boundary conditions (b, right hand side)
    b = assemble_periodic_boundary_conditions(p.boundary, panels, mesh)

    return InviscidSystem(A, b, mesh.panel_indices)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#

"""
    assemble_periodic_influence_matrix(v::Vortex, body_of_revolution, panels, mesh)

Assembles the "A" matrix (left hand side coefficient matrix).

# Arguments:
- `s::Singularity` : The singularity type used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

# Returns:
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_periodic_influence_matrix(::Singularity, panels, mesh) end

function assemble_periodic_influence_matrix(v::Vortex, panels, mesh)
    return assemble_periodic_vortex_matrix(v.order, panels, mesh)
end

"""
    assemble_periodic_vortex_matrix(::Order,  panels, mesh)

Assembles the coefficient matrix for a given order of singularity.

# Arguments:
- `o::Order` : The order of singularity used.
- `:Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

# Returns:
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_periodic_vortex_matrix(::Order, panels, mesh) end

function assemble_periodic_vortex_matrix(::Constant, panels, mesh)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies

    # initialize coefficient matrix
    TF = eltype(mesh.x)
    amat = zeros(TF, (N + nbodies, N + nbodies))

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panels --- ###
            for i in idx[m]
                for j in idx[n]

                    ### --- Calculate influence coefficient --- ###
                    s = j >= idx[m][end] - i ? -1.0 : 1.0
                    amat[i, j] =
                        s * calculate_periodic_vortex_influence(
                            Constant(), panels[m], panels[n], mesh, i, j
                        )
                end
            end

            if m == n

                ### --- Apply Back Substitution --- ###
                for i in idx[n]
                    sum = 0.0
                    jidx = idx[n][end] + 1 - i
                    for j in idx[m]
                        if j != jidx
                            sum += amat[j, i] * panels[n].panel_length[j]
                        end
                    end
                    dmagj = panels[n].panel_length[jidx]
                    amat[jidx, i] = -sum / dmagj
                end

                ### --- Apply Kutta Condition --- ###
                # put in the kutta condition for each airfoil (end rows of the system matrix)
                amat[N + m, idx[n][1]] = 1.0
                amat[N + m, idx[n][end]] = 1.0

                #put unit bound vortex value in each row
                amat[idx[m], idx[n][end] + n] .= 1.0
            end

            # # NOTE: this doesn't seem to change anything...
            # # Bound Vortex Correction
            # for i in 1:N
            #     for j in 1:M
            #         amat[i, j] += meshj.panels[j].length
            #     end
            # end

        end
    end

    return amat
end

#---------------------------------#
#    BOUNDARY CONDITION MATRIX    #
#---------------------------------#

"""
    assemble_periodic_boundary_conditions(::BoundaryCondition,  panels, mesh)

Assemble boundary condition vector.

# Arguments:
- `bc::BoundaryCondition` : The type of boundary condition to be used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

# Returns
 - `b::Matrix{Float}` : Boundary condition matrix
"""
function assemble_periodic_boundary_conditions(::BoundaryCondition, panels, mesh) end

function assemble_periodic_boundary_conditions(::Neumann, panels, mesh)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies

    # initialize boundary condition array
    TF = eltype(mesh.x)
    bc = zeros(TF, N + nbodies, 3)

    ### --- Loop through bodies --- ###
    for m in 1:nbodies

        # generate portion of boundary condition array associated with mth body
        bc[idx[m], 1] = [-cos(panels[m].panel_angle[i]) for i in idx[m]]
        bc[idx[m], 2] = [-sin(panels[m].panel_angle[i]) for i in idx[m]]
        bc[idx[m], 3] = [1.0 for i in idx[m]]
    end

    return bc
end
