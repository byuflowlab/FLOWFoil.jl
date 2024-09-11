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


