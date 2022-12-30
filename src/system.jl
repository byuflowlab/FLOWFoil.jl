#=

Inviscid System Functions

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                              GENERAL                               #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

abstract type System end

"""
    InviscidSystem

**Fields:**
 - `A::Array{Float,2}` : Coefficient Matrix on Left Hand Side.
 - `b::Array{Float,2}` : Boundary Condition Coefficient Vector on Right Hand Side.
 - `Ns::Array{Float}` : Array of numbers of nodes for each airfoil in the system.
"""
struct InviscidSystem{TA,TB,TI} <: System
    A::TA
    b::TB
    Ns::TI
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

"""
    generate_inviscid_system(problemtype::ProblemType, mesh, TEmesh)

**Arguments:**
- If PlanarProblem:
  - `problemtype::ProblemType` : ProblemType object for dispatch
  - `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
  - `mesh::Mesh` : Mesh for airfoil system to analyze.
  - `TEMesh::Mesh` : Trailing edge gap panel influence mesh.

- If AxisymmetricProblem:
  - `problemtype::ProblemType` : ProblemType object for dispatch
  - `body_of_revolution::Vector{Bool}` : flags whether bodies are bodies of revolution or not.
  - `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
  - `mesh::Mesh` : Mesh for airfoil system to analyze.

**Returns:**
` inviscid_system::InviscidSystem` : Inviscid System object containing influence and boundary condition matrices for the system.
"""
function generate_inviscid_system(problemtype::ProblemType, panels, mesh, TEmesh) end

######################################################################
#                                                                    #
#                               PLANAR                               #
#                                                                    #
######################################################################

#= NOTE:
The system assembly here is based on the Xfoil implementation, which does not strictly fit in the overall structure (it solves based on stream functions rather than potentials).
Likely, this will be moved elsewhere as other methods are developed.
=#
function generate_inviscid_system(pt::PlanarProblem, panels, mesh, TEmesh)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_influence_matrix(pt.singularity, mesh, TEmesh)

    # Get boundary conditions (b, right hand side)
    b = assemble_boundary_conditions(pt.boundary, panels, mesh, TEmesh)

    return InviscidSystem(A, b, mesh.node_indices)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#

"""
    assemble_influence_matrix(v::Vortex, mesh, TEmesh)

Assembles the "A" matrix (left hand side coefficient matrix).

**Arguments:**
- `s::Singularity` : The singularity type used.
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.
- `TEmesh::Mesh` : The mesh object associated with the trailing edge gap panels

**Returns:**
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_influence_matrix(::Singularity, mesh, TEmesh) end

function assemble_influence_matrix(v::Vortex, mesh, TEmesh)
    return assemble_vortex_matrix(v.order, mesh, TEmesh)
end

"""
    assemble_vortex_matrix(::Order, mesh, TEmesh)

Assembles the coefficient matrix for a given order of singularity.

**Arguments:**
- `o::Order` : The order of singularity used.
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.
- `TEmesh::Mesh` : The mesh object associated with the trailing edge gap panels

**Returns:**
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_vortex_matrix(::Order, mesh, TEmesh) end

function assemble_vortex_matrix(::Linear, mesh, TEmesh)

    # - Rename For Convenience - #
    pi = mesh.panel_indices
    ni = mesh.node_indices
    nb = mesh.nbodies

    # - Initialize Matrix - #
    TF = typeof(mesh.chord)
    amat = zeros(TF, ni[end][end] + nb, ni[end][end] + nb)

    ##### ----- Loop through the bodies being influenced ----- #####
    for m in 1:nb

        ##### ----- Loop through the bodies of influence ----- #####
        for n in 1:nb

            ### --- Populate main body of influence matrix --- ###
            for i in ni[m]
                for j in pi[n]
                    aij, aijp1 = calculate_vortex_influence(Linear(), mesh, i, j)

                    # add coefficients to matrix at correct nodes
                    if j == 1
                        amat[i, j + n - 1] = aij
                    else
                        amat[i, j + n - 1] += aij
                    end

                    amat[i, j + n] = aijp1
                end
            end

            ### --- Take Care of Matrix "Edges" --- ###
            #= Includes:
               - Kutta Condition
               - Trailing Edge Treatment
               - Constant Stream Function
            =#
            if n == m # when looking at the body influencing itself

                # - Find appropriate indices for this body - #

                # the ith row will be m after the nb row
                idx_i = ni[end][end] + m

                # the jth column for the first panel will be at the first index of the mth body
                idx_j1 = ni[m][1]

                # the jth column for the last panel will be at the last index of the mth body
                idx_j2 = ni[m][end]

                ### --- Apply Kutta Condition --- ###
                amat[idx_i, idx_j1] = 1.0
                amat[idx_i, idx_j2] = 1.0

                ### --- Insert Constant Stream Function Values --- ###
                #=NOTE:
                     We can re-use the indices we just found for the Kutta Condition,
                     but flip them so that they apply to the last columns instead.
                =#
                if !TEmesh.blunt_te[m]
                    #= If sharp trailing edge,
                      set the constant stream value
                      (associated with the Nth node equation for the airfoil)
                      to zero since that equation is replaced. =#
                    amat[idx_j1:(idx_j2 - 1), idx_i] .= -1.0

                    ### --- Trailing Edge Treatment --- ###
                    #= Replace last row of the submatrix with the extrapolation of the mean vortex strength to the trailing edge.
                    =#
                    # First zero out last row of submatrix
                    amat[ni[m][end], ni[m]] .= 0.0
                    # Then replace first and last elements in that row with extrapolation terms
                    amat[ni[m][end], ni[m][1]] = 1.0
                    amat[ni[m][end], ni[m][2]] = -2.0
                    amat[ni[m][end], ni[m][3]] = 1.0
                    amat[ni[m][end], ni[m][end] - 2] = -1.0
                    amat[ni[m][end], ni[m][end] - 1] = 2.0
                    amat[ni[m][end], ni[m][end]] = -1.0

                else
                    # otherwise keep everything at -1.0
                    amat[idx_j1:idx_j2, idx_i] .= -1.0

                    ### --- Add influence of trailing edge gap panel --- ###
                    for i in ni[m]

                        # Get panel influence coefficients
                        sigmate = calculate_source_influence(Constant(), TEmesh, i, m)
                        gammate = sum(calculate_vortex_influence(Linear(), TEmesh, i, m))

                        # Add/subtract from relevant matrix entries
                        amat[i, ni[m][1]] +=
                            0.5 * (gammate * TEmesh.tdp[m] - sigmate * TEmesh.txp[m])
                        amat[i, ni[m][end]] +=
                            0.5 * (sigmate * TEmesh.txp[m] - gammate * TEmesh.tdp[m])
                    end
                end
            end
        end
    end

    # println(size(amat))
    # display(amat[1:61, 1:61])
    # display(amat[62:122, 62:122])
    # display(amat[1:61, 62:122])
    # display(amat[62:122, 1:61])
    # display(amat[(end - 1):end, 1:61])
    # display(amat[(end - 1):end, 62:122])
    # display(amat[1:61, (end - 1):end])
    # display(amat[62:122, (end - 1):end])

    return amat
end

#---------------------------------#
#    BOUNDARY CONDITION MATRIX    #
#---------------------------------#

"""
    assemble_boundary_conditions(meshes)

Assemble boundary condition vector.

**Arguments:**
- `bc::BoundaryCondition` : The type of boundary condition to be used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.
- `TEmesh::Mesh` : The mesh object associated with the trailing edge gap panels

**Returns**
 - `b::Matrix{Float}` : Boundary condition matrix
"""
function assemble_boundary_conditions(::BoundaryCondition, panels, mesh, TEmesh) end

#= NOTE:
This implementation doesn't precisely fit.  As stated in other places, this Xfoil-like implementation will likely be moved with a better method can replace it.
=#
function assemble_boundary_conditions(::Dirichlet, panels, mesh, TEmesh)

    # - Rename For Convenience - #
    ni = mesh.node_indices
    nb = mesh.nbodies

    # initialize boundary condition array
    TF = typeof(mesh.chord)
    bmat = zeros(TF, ni[end][end] + nb, 2)

    ##### ----- Loop through system ----- #####
    for m in 1:nb

        ### --- Generate boundary condition array --- ###

        # if closed trailing edge, set second to last element of the boundary conditions to zero (for change in last node equation for sharp trailing edges)
        #=
          NOTE: mfoil does not do the following,
          but rather keeps the rhs as [-z,x] in all cases:
        =#
        bmat[ni[m][1]:ni[m][end - 1], 1] = [
            -mesh.nodes[i, 2] for i in ni[m][1]:ni[m][end - 1]
        ]

        bmat[ni[m][1]:ni[m][end - 1], 2] = [
            mesh.nodes[i, 1] for i in ni[m][1]:ni[m][end - 1]
        ]

        # if blunt trailing edge, no need for adjustment to last equation in submatrix.
        if TEmesh.blunt_te[m]
            bmat[ni[m][end], 1] = -mesh.nodes[ni[m][end], 2]
            bmat[ni[m][end], 2] = mesh.nodes[ni[m][end], 1]
        end
    end

    # display(bmat)

    return bmat
end

######################################################################
#                                                                    #
#                            AXISYMMETRIC                            #
#                                                                    #
######################################################################

function generate_inviscid_system(
    p::AxisymmetricProblem, panels::TP, mesh
) where {TP<:Panel}
    return generate_inviscid_system(p, [panels], mesh)
end

function generate_inviscid_system(p::AxisymmetricProblem, panels, mesh)
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

function assemble_ring_vortex_matrix(::Constant, body_of_revolution, panels, mesh)

    ### --- SETUP --- ###

    # Count number of bodies requiring a Kutta Condition
    nk = count(br -> br == false, body_of_revolution)
    kutta_count = 1

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies
    mesh2panel = mesh.mesh2panel

    # initialize coefficient matrix
    TF = eltype(mesh.m)
    amat = zeros(TF, (N + nk, N + nk))

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

            if m == n && !body_of_revolution[m]

                ### --- Apply Back Substitution --- ###
                for i in idx[n]
                    sum = 0.0
                    jidx = idx[n][end] + 1 - i
                    for j in idx[m]
                        if j != jidx
                            sum += amat[j, i] * panels[n].panel_length[mesh2panel[j]]
                        end
                    end
                    dmagj = panels[n].panel_length[mesh2panel[jidx]]
                    amat[jidx, i] = -sum / dmagj
                end

                ### --- Apply Kutta Condition --- ###
                # put in the kutta condition for each airfoil (end rows of the system matrix)
                amat[N + kutta_count, idx[n][1]] = 1.0
                amat[N + kutta_count, idx[n][end]] = 1.0

                #put unit bound vortex value in each row
                amat[idx[m], idx[n][end] + kutta_count] .= 1.0

                kutta_count += 1
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

function assemble_ring_boundary_conditions(::Neumann, body_of_revolution, panels, mesh)

    ### --- SETUP --- ###

    # Count number of bodies requiring a Kutta Condition
    nk = count(br -> br == false, body_of_revolution)

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies
    mesh2panel = mesh.mesh2panel

    # initialize boundary condition array
    TF = eltype(mesh.m)
    bc = zeros(TF, N + nk)

    ### --- Loop through bodies --- ###
    for m in 1:nbodies

        # generate portion of boundary condition array associated with mth body
        if body_of_revolution[m]
            bc[idx[m], 1] = [-cos(panels[m].panel_angle[mesh2panel[i]]) for i in idx[m]]
        else
            bc[idx[m], 1] = [
                1.0 - cos(panels[m].panel_angle[mesh2panel[i]]) for i in idx[m]
            ]
        end
    end

    return bc
end

######################################################################
#                                                                    #
#                             PERIODIC                               #
#                                                                    #
######################################################################

function generate_inviscid_system(p::PeriodicProblem, panels::TP, mesh) where {TP<:Panel}
    return generate_inviscid_system(p, [panels], mesh)
end

function generate_inviscid_system(p::PeriodicProblem, panels, mesh)
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

**Arguments:**
- `s::Singularity` : The singularity type used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

**Returns:**
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_periodic_influence_matrix(::Singularity, panels, mesh) end

function assemble_periodic_influence_matrix(v::Vortex, panels, mesh)
    return assemble_periodic_vortex_matrix(v.order, panels, mesh)
end

"""
    assemble_periodic_vortex_matrix(::Order,  panels, mesh)

Assembles the coefficient matrix for a given order of singularity.

**Arguments:**
- `o::Order` : The order of singularity used.
- `:Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

**Returns:**
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

**Arguments:**
- `bc::BoundaryCondition` : The type of boundary condition to be used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

**Returns**
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
