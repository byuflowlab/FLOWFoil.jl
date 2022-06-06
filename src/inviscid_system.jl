#=
Inviscid System Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    assemble_vortex_coefficients(meshi, meshj)

Assemble matrix of vortex strength coefficients.

This function only assembles the NxM portion of the system influence coefficient matrix associated with the M-1 panels of meshj acting on the N nodes of meshi. It does not include the kutta condition or the influence of the constant stream function on the airfoil nodes.

**Arguments:**
 - 'meshi::BodyMesh' : mesh being influenced.
 - 'meshj::BodyMesh' : mesh doing the influencing.

"""
function assemble_vortex_coefficients(mehsi, meshj)

    # get nodes for convenience
    nodesi = meshi.airfoil_nodes
    nodesj = meshj.airfoil_nodes

    # get system size
    N = length(nodesi)
    M = length(nodesj)

    #initialize NxN coefficient matrix
    amat = [0.0 for i in 1:N, j in 1:M]

    # loop through setting up influence coefficients
    for i in 1:N
        for j in 1:(M - 1)

            # obtain influence coefficient for ith evaluation point and j and j+1 panel
            aij, aijp1 = FLOWFoil.get_vortex_influence(nodesj[j], nodesj[j + 1], nodesi[i])

            # add coefficients to matrix at correct nodes
            if j == 1
                amat[i, j] = aij
            else
                amat[i, j] += aij
            end

            amat[i, j + 1] = aijp1
        end
    end

    # add in trailing edge contributions
    # NOTE!: mfoil applies these all the time.
    if meshj.blunt_te
        for i in 1:N

            # Get panel influence coefficients
            sigmate = get_source_influence(nodesj[M], nodesj[1], nodesi[i])
            gammate = sum(get_vortex_influence(nodesj[M], nodesj[1], nodesi[i]))

            # Add/subtract from relevant matrix entries
            amat[i, 1] += 0.5 * (gammate * mesh.tdp - sigmate * mesh.txp)
            amat[i, M] += 0.5 * (sigmate * mesh.txp - gammate * mesh.tdp)
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

    return amat
end

"""
    assemble_vortex_matrix(meshes)

Assemble vortex coefficient matrix with full N+n x N+n system, including Kutta condition, where n is the number of meshes (airfoils) in the system, and N is the total number of nodes between all the meshes.

**Arguments:**
 - 'meshes::Array{BodyMesh}' : Mesh System for which to solve.

"""
function assemble_vortex_matrix(meshes)

    # size sysetm
    N, Ns = size_system(meshes)
    n = length(Ns)
    offset = get_offset(Ns)

    # initialize coefficient matrix
    amat = [0.0 for i in 1:(N + n), j in 1:(N + n)]

    # Loop through system
    for x in 1:n
        for y in 1:n

            # get influence coefficients for each portion of the sysetm (mesh y acts on mesh x)
            # put things in the correct place in the system matrix
            amat[(1 + offset[x]):(Ns[x] + offset[x]), (1 + offset[y]):(Ns[y] + offset[y])] .= assemble_vortex_coefficients(
                meshes[x], meshes[y]
            )

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
            end
        end
    end

    return amat
end

"""
    assemble_boundary_conditions(meshes)

Assemble boundary condition vector.

**Arguments:**
 - 'meshes::Array{BodyMesh}' : mesh system for which to solve

**Returns**
 - psi_inf::Array{Float,2} : Boundary condition array.
"""
function assemble_boundary_conditions(meshes)

    # size the system
    N, Ns = size_system(meshes)
    offset = get_offset(Ns)

    # initialize boundary condition array
    bc = [0.0 for i in 1:N, j in 1:2]

    # Loop through system
    for m in 1:length(Ns)

        # get node locations for convenience
        nodes = mesh.airfoil_nodes
        N = length(nodes)

        # generate boundary condition array
        if !mesh.blunt_te
            # NOTE: mfoil does not do the following, but rather keeps the rhs as [-z,x] in all cases:
            # if closed trailing edge, set last element of the boundary conditions to zero
            bc[(1 + offset[m]):(Ns[m] + offset[m] - 1), :] .= [
                [-nodes[i][2]; nodes[i][1]] for i in 1:N
            ]

        else
            bc[(1 + offset[m]):(Ns[m] + offset[m]), :] .= [
                [-nodes[i][2]; nodes[i][1]] for i in 1:N
            ]
        end
    end

    return bc
end

"""
    get_inviscid_system(mesh)

Calculate, then gather the vortex and boundary condition matricies into an InviscidSystem object.

**Arguments:**
- 'mesh::BodyMesh' : BodyMesh for airfoil to analyze.
"""
function get_inviscid_system(mesh)

    # Get coeffiecient matrix (A, left hand side)
    vcoeffmat = assemble_vortex_matrix(mesh)

    # Get boundary conditions (RHS)
    bccoeffvec = assemble_boundary_conditions(mesh)

    return InviscidSystem(vcoeffmat, bccoeffvec)
end
