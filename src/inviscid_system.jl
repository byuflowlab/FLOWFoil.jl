#=
Inviscid System Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    assemble_vortex_coefficients(mesh)

Assemble matrix of vortex strength coefficients.

This function only assembles the NxN portion of the influence coefficient matrix. It does not include the kutta condition.
Also note that multibody capabilities have not yet been implemented.

**Arguments:**
 - 'mesh::BodyMesh' : mesh system for which to find influence coefficient matrix.

"""
function assemble_vortex_coefficients(mesh)

    # get nodes for convenience
    nodes = mesh.airfoil_nodes

    # get system size
    N = length(nodes)

    #initialize NxN coefficient matrix
    amat = Array{Float64,2}(undef, N, N)

    # loop through setting up influence coefficients
    for i in 1:N
        for j in 1:(N - 1)

            # obtain influence coefficient for ith evaluation point and j and j+1 panel
            aij, aijp1 = FLOWFoil.get_vortex_influence(nodes[j], nodes[j + 1], nodes[i])

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
    if mesh.blunt_te
        for i in 1:N

            # Get panel influence coefficients
            sigmate = get_source_influence(nodes[N], nodes[1], nodes[i])
            gammate = sum(get_vortex_influence(nodes[N], nodes[1], nodes[i]))
            # println("i: $i")

            # Add/subtract from relevant matrix entries
            amat[i, 1] += 0.5 * (gammate * mesh.tdp - sigmate * mesh.txp)
            amat[i, N] += 0.5 * (sigmate * mesh.txp - gammate * mesh.tdp)
        end
    end

    if !mesh.blunt_te
        # Replace Nth row of the matrix with the extrapolation of the mean vortex strength to the trailing edge.
        amat[end, :] .= 0.0
        amat[end, 1] = 1.0
        amat[end, 2] = -2.0
        amat[end, 3] = 1.0
        amat[end, end - 2] = -1.0
        amat[end, end - 1] = 2.0
        amat[end, end] = -1.0
    end

    return amat
end

"""
    assemble_vortex_matrix(mesh)

Assemble vortex coefficient matrix with full N+1 x N+1 system, including Kutta condition.

**Arguments:**
 - 'mesh::BodyMesh' : Mesh System for which to solve.

"""
function assemble_vortex_matrix(mesh)

    #get NxN square of a matrix
    amat = assemble_vortex_coefficients(mesh)

    # get size of amat:
    r, c = size(amat)

    # Add column of ones to end of matrix for Psi_0
    psi0 = -1.0 * ones(r)
    # if sharp trailing edge, set a[N,N+1] to zero
    if !mesh.blunt_te
        psi0[end] = 0.0
    end
    amat = [amat psi0]

    # add kutta condition row to bottom of matrix
    kutta = zeros(c + 1)
    kutta[1] = 1.0
    kutta[end - 1] = 1.0
    amat = [amat; kutta']

    return amat
end

"""
    assemble_boundary_conditions(mesh, freestream)

Assemble boundary condition vector.

**Arguments:**
 - 'mesh::BodyMesh' : mesh system for which to solve

**Returns**
 - output::type : description.
"""
function assemble_boundary_conditions(mesh)

    # get node locations for convenience
    nodes = mesh.airfoil_nodes
    N = length(nodes)

    # generate boundary condition array
    psi_inf1 = [-nodes[i][2] for i in 1:N]
    psi_inf2 = [nodes[i][1] for i in 1:N]
    psi_inf = [psi_inf1 psi_inf2]

    #NOTE: mfoil does not do the following, but rather keeps the rhs as [-z,x] in all cases:
    # if closed trailing edge, set last element of psi_inf to zero
    if !mesh.blunt_te
        psi_inf[end, :] = [0.0 0.0]
    end

    # For Kutta Condition, RHS = 0.0
    psi_inf = vcat(psi_inf, [0.0 0.0])

    return psi_inf
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
