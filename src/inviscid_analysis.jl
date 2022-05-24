#=
Inviscid System Functions

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

#TODO: probably want to create another version of this that takes in a mesh only.
"""
    assemblevortexcoefficients(meshsystem)

Assemble matrix of vortex strength coefficients.

This function only assembles the NxN portion of the influence coefficient matrix. It does not include the kutta condition.
Also note that multibody capabilities have not yet been implemented.

**Arguments:**
 - 'meshsystem::MeshSystem' : mesh system for which to find influence coefficient matrix.

"""
function assemblevortexcoefficients(meshsystem)

    # get system size
    N, Ns = FLOWFoil.sizesystem(meshsystem)

    #check that multibody system hasn't been used yet
    if N != Ns[1]
        @warn("Multi-body systems are not yet supported")
    end

    # set N to Ns[1] for now.
    N = Ns[1]

    # get nodes for convenience
    nodes = meshsystem.meshes[1].airfoil_nodes

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
    #!NOTE: mfoil seems to apply blunt trailing edge contributions always, whether or not there actually is a blunt trailing edge.  Theoretically, the TE contributions would go to zero if the TE points were coincident, but they may, in fact, be within 1e-10.
    # TODO: Probably move most of these calculations to a separate function for organizational purposes.
    # get bisection vector
    # get vector along first panel
    if meshsystem.meshes[1].blunt_te
        nd1, _ = get_d(nodes[2], nodes[1])

        # get vector along second panel
        dn, _ = get_d(nodes[end - 1], nodes[end])

        # calculate vector that bisects the first and last panel vectors using formula c = |a|*b + |b|*a
        bisector = nd1 * sqrt(dn[1]^2 + dn[2]^2) + dn * sqrt(nd1[1]^2 + nd1[2]^2)

        # normalize to get the unit vector
        ttehat = bisector / sqrt(bisector[1]^2 + bisector[2]^2)

        # get panel vector
        dte, _ = get_d(nodes[end], nodes[1])

        # normalize panelvector
        dtehat = dte / sqrt(dte[1]^2 + dte[2]^2)

        # get dot product of bisection vector and panel vector.
        tdp = ttehat[1] * dtehat[1] + ttehat[2] * dtehat[2]

        # get cross product of bisection vector and panel vector
        txp = abs(ttehat[1] * dtehat[2] - ttehat[2] * dtehat[1])

        # Apply trailing edge contributions
        for i in 1:N

            # Get panel influence coefficients
            sigmate = get_source_influence(nodes[N], nodes[1], nodes[i])
            gammate = sum(get_vortex_influence(nodes[N], nodes[1], nodes[i]))

            # Add/subtract from relevant matrix entries
            amat[i, 1] += 0.5 * (gammate * tdp - sigmate * txp)
            amat[i, N] += 0.5 * (sigmate * txp - gammate * tdp)
        end
    else
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
    assemblematrixa(meshsystem)

Assemble vortex coefficient matrix with full N+1 x N+1 system, including Kutta condition.

**Arguments:**
 - 'meshsystem:MeshSystem' : Mesh System for which to solve.

"""
function assemblematrixa(meshsystem)

    #get NxN square of a matrix
    amat = assemblevortexcoefficients(meshsystem)

    # get size of amat:
    r, c = size(amat)

    # Add column of ones to end of matrix for Psi_0
    psi0 = -1.0 * ones(r)
    # if sharp trailing edge, set a[N,N+1] to zero
    if !meshsystem.meshes[1].blunt_te
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
    assembleboundaryconditions(meshsystem, freestream)

Assemble boundary condition vector.

Note that multiple operation conditions is not yet supported.
If freestream fields contain more than one item, only the first will be used.

**Arguments:**
 - 'meshsystem::MeshSystem' : mesh system for which to solve
 - 'freestream::Freestream' : freestream parameters

**Returns**
 - output::type : description.
"""
function assembleboundaryconditions(meshsystem)

    # get node locations for convenience
    nodes = meshsystem.meshes[1].airfoil_nodes
    N = length(nodes)

    # generate boundary condition array
    psi_inf1 = [-nodes[i][2] for i in 1:N]
    psi_inf2 = [nodes[i][1] for i in 1:N]
    psi_inf = [psi_inf1 psi_inf2]

    #NOTE: mfoil does not do the following, but rather keeps the rhs as [-z,x] in all cases:
    # if closed trailing edge, set last element of psi_inf to zero
    if !meshsystem.meshes[1].blunt_te
        psi_inf[end, :] = [0.0 0.0]
    end

    # For Kutta Condition, RHS = 0.0
    psi_inf = vcat(psi_inf, [0.0 0.0])

    return psi_inf
end
