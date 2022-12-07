#=

Inviscid System Functions

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                           GENERAL TYPES                            #
#                                                                    #
######################################################################

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

######################################################################
#                                                                    #
#                    SYSTEM GENERATION FUNCTIONS                     #
#                                                                    #
######################################################################

"""
    generate_inviscid_system(problemtype::ProblemType, mesh, TEmesh)

**Arguments:**
- `problemtype::ProblemType` : ProblemType object for dispatch
- `mesh::Mesh` : Mesh for airfoil system to analyze.
- `TEMesh::Mesh` : Trailing edge gap panel influence mesh.

**Returns:**
` inviscid_system::InviscidSystem` : Inviscid System object containing influence and boundary condition matrices for the system.
"""
function generate_inviscid_system(problemtype::ProblemType, panels, mesh, TEmesh)
    return generate_inviscid_system(problemtype.method, panels, mesh, TEmesh)
end

#---------------------------------#
#             PLANAR              #
#---------------------------------#

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

"""
    assemble_influence_matrix(v::Vortex, mesh, TEmesh)

Assembles the "A" matrix (left hand side coefficient matrix.

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
                        amat[i, j] = aij
                    else
                        amat[i, j] += aij
                    end

                    amat[i, j + 1] = aijp1
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

                else
                    # otherwise keep everything at -1.0
                    amat[idx_j1:idx_j2, idx_i] .= -1.0

                    ### --- Add influence of trailing edge gap panel --- ###
                    for i in ni[m]

                        # Get panel influence coefficients
                        sigmate = get_source_influence(Constant(), TEmesh, i, m)
                        gammate = sum(get_vortex_influence(Constant(), TEmesh, i, m))

                        # Add/subtract from relevant matrix entries
                        amat[i, 1] +=
                            0.5 * (gammate * TEmesh.tdp[m] - sigmate * TEmesh.txp[m])
                        amat[i, M] +=
                            0.5 * (sigmate * TEmesh.txp[m] - gammate * TEmesh.tdp[m])
                    end
                end

                ### --- Trailing Edge Treatment --- ###
                #= Replace last row of the submatrix with the extrapolation of the mean vortex strength to the trailing edge.
                =#
                # First zero out last row of submatrix
                amat[ni[m][end], :] .= 0.0
                # Then replace first and last elements in that row with extrapolation terms
                amat[ni[m][end], 1] = 1.0
                amat[ni[m][end], 2] = -2.0
                amat[ni[m][end], 3] = 1.0
                amat[ni[m][end], ni[n][end] - 2] = -1.0
                amat[ni[m][end], ni[n][end] - 1] = 2.0
                amat[ni[m][end], ni[n][end]] = -1.0
            end
        end
    end

    return amat
end

"""
    assemble_boundary_conditions(meshes)

Assemble boundary condition vector.

**Arguments:**
- `bc::BoundaryCondition` : The type of boundary condition to be used.
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
            -panels.panel_edges[i, 1, 2] for i in ni[m][1]:ni[m][end - 1]
        ]
        bmat[ni[m][1]:ni[m][end - 1], 2] = [
            panels.panel_edges[i, 1, 1] for i in ni[m][1]:ni[m][end - 1]
        ]

        # if blunt trailing edge, no need for adjustment to last equation in submatrix.
        if TEmesh.blunt_te[m]
            bmat[ni[m][end], 1] = -panels.panel_edges[i, 2, 2]
            bmat[ni[m][end], 2] = panels.panel_edges[i, 2, 1]
        end
    end

    return bmat
end

#---------------------------------#
#           AXISYMMETRIC          #
#---------------------------------#

"""
**Arguments:**
- `mesh::Array{AxisymMesh}` : AxisymMesh for airfoil to analyze.
"""
function generate_inviscid_system(::AxisymmetricProblem, mesh)
    # Get coeffiecient matrix (A, left hand side)
    A, Ns = assemble_vortex_matrix(mesh)

    # Get boundary conditions (b, right hand side)
    b = assemble_boundary_conditions(mesh)

    return InviscidSystem(A, b, Ns)
end
