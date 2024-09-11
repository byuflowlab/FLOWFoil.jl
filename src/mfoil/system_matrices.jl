function generate_system_matrices(pt::Mfoil, panels, mesh, TEmesh)
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
    pidx = mesh.panel_indices
    nidx = mesh.node_indices
    nb = mesh.nbodies

    # - Initialize Matrix - #
    TF = typeof(mesh.chord)
    amat = zeros(TF, nidx[end][end] + nb, nidx[end][end] + nb)

    ##### ----- Loop through the bodies being influenced ----- #####
    for m in 1:nb

        ##### ----- Loop through the bodies of influence ----- #####
        for n in 1:nb

            ### --- Populate main body of influence matrix --- ###
            for i in nidx[m]
                for j in pidx[n]
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
                idx_i = nidx[end][end] + m

                # the jth column for the first panel will be at the first index of the mth body
                idx_j1 = nidx[m][1]

                # the jth column for the last panel will be at the last index of the mth body
                idx_j2 = nidx[m][end]

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
                    amat[nidx[m][end], nidx[m]] .= 0.0
                    # Then replace first and last elements in that row with extrapolation terms
                    amat[nidx[m][end], nidx[m][1]] = 1.0
                    amat[nidx[m][end], nidx[m][2]] = -2.0
                    amat[nidx[m][end], nidx[m][3]] = 1.0
                    amat[nidx[m][end], nidx[m][end] - 2] = -1.0
                    amat[nidx[m][end], nidx[m][end] - 1] = 2.0
                    amat[nidx[m][end], nidx[m][end]] = -1.0

                else
                    # otherwise keep everything at -1.0
                    amat[idx_j1:idx_j2, idx_i] .= -1.0

                    ### --- Add influence of trailing edge gap panel --- ###
                    for i in nidx[m]

                        # Get panel influence coefficients
                        sigmate = calculate_source_influence(Constant(), TEmesh, i, m)
                        gammate = sum(calculate_vortex_influence(Linear(), TEmesh, i, m))

                        # Add/subtract from relevant matrix entries
                        amat[i, nidx[m][1]] +=
                            0.5 * (gammate * TEmesh.tdp[m] - sigmate * TEmesh.txp[m])
                        amat[i, nidx[m][end]] +=
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
    nidx = mesh.node_indices
    nb = mesh.nbodies

    # initialize boundary condition array
    TF = typeof(mesh.chord)
    bmat = zeros(TF, nidx[end][end] + nb, 2)

    ##### ----- Loop through system ----- #####
    for m in 1:nb

        ### --- Generate boundary condition array --- ###

        # if closed trailing edge, set second to last element of the boundary conditions to zero (for change in last node equation for sharp trailing edges)
        #=
          NOTE: mfoil does not do the following,
          but rather keeps the rhs as [-z,x] in all cases:
        =#
        bmat[nidx[m][1]:nidx[m][end - 1], 1] = [
            -mesh.nodes[i, 2] for i in nidx[m][1]:nidx[m][end - 1]
        ]

        bmat[nidx[m][1]:nidx[m][end - 1], 2] = [
            mesh.nodes[i, 1] for i in nidx[m][1]:nidx[m][end - 1]
        ]

        # if blunt trailing edge, no need for adjustment to last equation in submatrix.
        if TEmesh.blunt_te[m]
            bmat[nidx[m][end], 1] = -mesh.nodes[nidx[m][end], 2]
            bmat[nidx[m][end], 2] = mesh.nodes[nidx[m][end], 1]
        end
    end

    # display(bmat)

    return bmat
end
