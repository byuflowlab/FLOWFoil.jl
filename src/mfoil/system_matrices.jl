#=
NOTE: THERE ARE 2 IMPLEMENTATIONS HERE, NOT SURE WHICH ONE IS THE RIGHT ONE AT THIS POINT
probably the second set in this file (after major header)
=#

"""
    assemble_vortex_coefficients(meshi, meshj)

Assemble matrix of vortex strength coefficients.

This function only assembles the NxM portion of the system influence coefficient matrix associated with the M-1 panels of meshj acting on the N nodes of meshi. It does not include the kutta condition or the influence of the constant stream function on the airfoil nodes.

**Arguments:**
 - `meshi::PlanarMesh` : mesh being influenced.
 - `meshj::PlanarMesh` : mesh doing the influencing.
 - `trailing_edge_treatment::Bool` : flag for whether to treat trailing edge or not (is meshi==meshj?)
"""
function assemble_vortex_coefficients(meshi, meshj, args...)

    # get nodes for convenience
    nodesi = meshi.nodes
    nodesj = meshj.nodes

    # get system size
    N = length(nodesi)
    M = length(nodesj)

    rtype = promote_type(eltype(nodesi[1]), eltype(nodesj[1]))

    amat = zeros(rtype, N, M)

    return assemble_vortex_coefficients!(amat, meshi, meshj, args...)
end
function assemble_vortex_coefficients!(amat, meshi, meshj, trailing_edge_treatment)

    # get nodes for convenience
    nodesi = meshi.nodes
    nodesj = meshj.nodes

    # get system size
    N = length(nodesi)
    M = length(nodesj)

    #initialize NxN coefficient matrix
    # amat = [0.0 for i in 1:N, j in 1:M]

    # loop through setting up influence coefficients
    for i in 1:N
        for j in 1:(M - 1)

            # obtain influence coefficient for ith evaluation point and j and j+1 panel
            # aij, aijp1 = get_vortex_influence(nodesj[j], nodesj[j + 1], nodesi[i])

            # NOTE: Here we add a little offset to keep ForwardDiff from
            #   returning NaN on the self-influence (case sqrt(0) when
            #   the target is one of the nodes)
            aij, aijp1 = get_vortex_influence(nodesj[j], nodesj[j + 1], nodesi[i] .+ 1e-9)

            # add coefficients to matrix at correct nodes
            if j == 1
                amat[i, j] = aij
            else
                amat[i, j] += aij
            end

            amat[i, j + 1] = aijp1
        end
    end

    if trailing_edge_treatment
        # add in trailing edge contributions
        # NOTE!: mfoil applies these all the time.
        if meshj.blunt_te
            for i in 1:N

                # Get panel influence coefficients
                sigmate = get_source_influence(nodesj[M], nodesj[1], nodesi[i])

                # gammate = sum(get_vortex_influence(nodesj[M], nodesj[1], nodesi[i]))

                # NOTE: Here we add a little offset to keep ForwardDiff from
                #   returning NaN on the self-influence (case sqrt(0) when
                #   the target is one of the nodes)
                gammate = sum(get_vortex_influence(nodesj[M], nodesj[1], nodesi[i] .+ 1e-9))

                # Add/subtract from relevant matrix entries
                amat[i, 1] += 0.5 * (gammate * meshj.tdp - sigmate * meshj.txp)
                amat[i, M] += 0.5 * (sigmate * meshj.txp - gammate * meshj.tdp)
            end

        else
            # Replace Nth row of the matrix with the extrapolation of the mean vortex strength to the trailing edge.
            amat[N, :] .= 0.0
            amat[N, 1] = 1.0
            amat[N, 2] = -2.0
            amat[N, 3] = 1.0
            amat[N, M - 2] = -1.0
            amat[N, M - 1] = 2.0
            amat[N, M] = -1.0
        end
    end

    return amat
end

"""
    assemble_vortex_matrix(meshes)

Assemble vortex coefficient matrix with full N+n x N+n system, including Kutta condition, where n is the number of meshes (airfoils) in the system, and N is the total number of nodes between all the meshes.

**Arguments:**
 - `meshes::Array{PlanarMesh}` : Mesh System for which to solve.
"""
function assemble_vortex_matrix(meshes)

    # size sysetm
    N, Ns = size_system(meshes)
    n = length(Ns)

    rtype = promote_type(eltype.(mesh.nodes[1] for mesh in meshes)...)

    # initialize coefficient matrix
    amat = zeros(rtype, N + n, N + n)

    return assemble_vortex_matrix!(amat, meshes)
end
function assemble_vortex_matrix!(amat, meshes)

    # size sysetm
    N, Ns = size_system(meshes)
    n = length(Ns)
    offset = get_offset(Ns)

    # # initialize coefficient matrix
    # amat = [0.0 for i in 1:(N + n), j in 1:(N + n)]

    # Loop through system
    for x in 1:n
        for y in 1:n
            if x == y
                # put in the kutta condition for each airfoil (end rows of the system matrix)
                amat[N + x, 1 + offset[y]] = 1.0
                amat[N + x, Ns[y] + offset[y]] = 1.0

                # put in the constant stream function influences for each airfoil (end columns of the system matrix)
                if !meshes[x].blunt_te
                    # if sharp trailing edge, set the constant stream value associated with the Nth node equation for the airfoil to zero (since that equation is replaced)
                    amat[(1 + offset[x]):(Ns[x] + offset[x] - 1), N + y] .= -1.0
                else
                    # otherwise keep everything at -1.0
                    amat[(1 + offset[x]):(Ns[x] + offset[x]), N + y] .= -1.0
                end
                trailing_edge_treatment = true
            else
                trailing_edge_treatment = false
            end

            # get influence coefficients for each portion of the sysetm (mesh y acts on mesh x)
            # put things in the correct place in the system matrix
            amat[(1 + offset[x]):(Ns[x] + offset[x]), (1 + offset[y]):(Ns[y] + offset[y])] .= assemble_vortex_coefficients(
                meshes[x], meshes[y], trailing_edge_treatment
            )
        end
    end

    return amat, Ns
end

"""
    assemble_boundary_conditions(meshes)

Assemble boundary condition vector.

**Arguments:**
 - `meshes::Array{PlanarMesh}` : mesh system for which to solve

**Returns**
 - `psi_inf::Array{Float,2}` : Boundary condition array.
"""
function assemble_boundary_conditions(meshes)

    # size the system
    N, Ns = size_system(meshes)

    rtype = promote_type(eltype.(mesh.nodes[1] for mesh in meshes)...)

    # initialize boundary condition array
    bc = zeros(rtype, N + length(Ns), 1:2)

    return assemble_boundary_conditions!(bc, meshes)
end

function assemble_boundary_conditions!(bc, meshes)

    # size the system
    N, Ns = size_system(meshes)
    offset = get_offset(Ns)

    # # initialize boundary condition array
    # bc = [0.0 for i in 1:(N + length(Ns)), j in 1:2]

    # Loop through system
    for m in 1:length(Ns)

        # get node locations for convenience
        nodes = meshes[m].nodes
        N = length(nodes)

        # generate boundary condition array
        if !meshes[m].blunt_te
            # NOTE: mfoil does not do the following, but rather keeps the rhs as [-z,x] in all cases:
            # if closed trailing edge, set last element of the boundary conditions to zero
            bc[(1 + offset[m]):(Ns[m] + offset[m] - 1), 1] = [
                -nodes[i][2] for i in 1:(N - 1)
            ]
            bc[(1 + offset[m]):(Ns[m] + offset[m] - 1), 2] = [
                nodes[i][1] for i in 1:(N - 1)
            ]

        else
            bc[(1 + offset[m]):(Ns[m] + offset[m]), 1] = [-nodes[i][2] for i in 1:N]
            bc[(1 + offset[m]):(Ns[m] + offset[m]), 2] = [nodes[i][1] for i in 1:N]
        end
    end

    return bc
end

######################################################################
#                                                                    #
#                        OTHER IMPLEMENTATION                        #
#                                                                    #
######################################################################
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

