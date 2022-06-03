#=
Inviscid System Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    assemble_vortex_coefficients(meshes)

Assemble matrix of vortex strength coefficients.

This function only assembles the NxN portion of the influence coefficient matrix. It does not include the kutta condition.
Also note that multibody capabilities have not yet been implemented.

**Arguments:**
 - 'meshes::BodyMesh' : mesh system for which to find influence coefficient matrix.

"""
function assemble_vortex_matrix(meshes)

    # size system
    N, Ns = size_system(meshes)
    M = length(meshes)

    #initialize NxN coefficient matrix
    amat = [0.0 for i in 1:(N + M), j in 1:(N + M)]

    # Constant Stream Function Component
    amat[1:(end - M), (N + 1):end] .= -1.0

    # set offset values
    offset = [0; cumsum(Ns[1:(end - 1)])]

    # loop through setting up influence coefficients
    # for each airfoil
    for m in 1:M

        # get nodes for convenience
        p = meshes[m].airfoil_nodes

        #for each airfoil
        for n in 1:M

            # get nodes for convenience
            q = meshes[n].airfoil_nodes

            # for each node in airfoil m
            for i in 1:Ns[m]

                # for each panel in airfoil n
                for j in 1:(Ns[n] - 1)

                    # obtain influence coefficient for ith evaluation point and j and j+1 panel
                    aij, aijp1 = FLOWFoil.get_vortex_influence(q[j], q[j + 1], p[i])

                    # add coefficients to matrix at correct nodes
                    amat[i + offset[m], j + offset[n]] += aij
                    amat[i + offset[m], j + 1 + offset[n]] = aijp1
                end #for panels

                # if airfoil m has a blunt trailing edge, add in TE panel influences
                if meshes[n].blunt_te
                    # Get panel influence coefficients
                    sigmate = get_source_influence(q[end], q[1], p[i])
                    gammate = sum(get_vortex_influence(q[end], q[1], p[i]))

                    # Add/subtract from relevant matrix entries
                    amat[i + offset[m], 1 + offset[n]] +=
                        0.5 * (gammate * meshes[n].tdp - sigmate * meshes[n].txp)
                    amat[i + offset[m], Ns[n] + offset[n]] +=
                        0.5 * (sigmate * meshes[n].txp - gammate * meshes[n].tdp)
                end
            end # for nodes

            # Kutta Condition
            amat[end - (M - m), 1 + offset[n]] = 1.0
            amat[end - (M - m), Ns[n] + offset[n]] = 1.0

            # if airfoil m has a sharp trailing edge, and m==n, swap out last line
            if !meshes[m].blunt_te && m == n

                # set constant stream value to zero on last panel
                amat[Ns[m] + offset[m], end - (M - n)] = 0.0
                # Replace Nth row of the matrix with the extrapolation of the mean vortex strength to the trailing edge.
                amat[Ns[m] + offset[m], (1 + offset[n]):(Ns[n] + offset[n])] .= 0.0
                amat[Ns[m] + offset[m], 1 + offset[n]] = 1.0
                amat[Ns[m] + offset[m], 2 + offset[n]] = -2.0
                amat[Ns[m] + offset[m], 3 + offset[n]] = 1.0
                amat[Ns[m] + offset[m], Ns[n] + offset[n] - 2] = -1.0
                amat[Ns[m] + offset[m], Ns[n] + offset[n] - 1] = 2.0
                amat[Ns[m] + offset[m], Ns[n] + offset[n]] = -1.0
            end #if blunt TE
        end # for airfoil n
    end # for airfoil n

    return amat, Ns
end

"""
    assemble_boundary_conditions(mesh, freestream)

Assemble boundary condition vector.

**Arguments:**
 - 'mesh::BodyMesh' : mesh system for which to solve

**Returns**
 - output::type : description.
"""
function assemble_boundary_conditions(meshes)

    # size system
    N, Ns = size_system(meshes)
    offset = [0; cumsum(Ns[1:(end - 1)])]
    M = length(meshes)

    # Initialize Boundary Condition Matrix
    psi_inf = [0.0 for i in 1:(N + M), j in 1:2]

    for m in 1:M
        # Populate boundary conditions
        psi_inf[(1 + offset[m]):(Ns[m] + offset[m]), 1] = [
            -meshes[m].airfoil_nodes[i][2] for i in 1:Ns[m]
        ]
        psi_inf[(1 + offset[m]):(Ns[m] + offset[m]), 2] = [
            meshes[m].airfoil_nodes[i][1] for i in 1:Ns[m]
        ]

        # NOTE: mfoil does not do the following, but rather keeps the rhs as [-z,x] in all cases:
        # if closed trailing edge, set last element of psi_inf to zero
        if !meshes[m].blunt_te
            psi_inf[Ns[m] + offset[m], :] = [0.0 0.0]
        end
    end

    return psi_inf
end

"""
    get_inviscid_system(mesh)

Calculate, then gather the vortex and boundary condition matricies into an InviscidSystem object.

**Arguments:**
- 'mesh::BodyMesh' : BodyMesh for airfoil to analyze.
"""
function get_inviscid_system(meshes)

    # Get coeffiecient matrix (A, left hand side)
    vcoeffmat, Ns = assemble_vortex_matrix(meshes)

    # Get boundary conditions (RHS)
    bccoeffvec = assemble_boundary_conditions(meshes)

    return InviscidSystem(vcoeffmat, bccoeffvec, Ns)
end
