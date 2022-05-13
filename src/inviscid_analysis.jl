#=
Inviscid Panel Method Solver

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\overbar{\\Psi}^\\gamma\$

**Arguments:**
 - 'theta1::Float' : angle between panel and vector from node1 to evaluation point
 - 'theta2::Float' : angle between panel and vector from node2 to evaluation point
 - 'ln1::Float' : value of ln(rmag1), which may be that or 0.0, depening on evaluation point location
 - 'ln2::Float' : value of ln(rmag2), which may be that or 0.0, depening on evaluation point location
 - 'dmag::Float' : panel length
 - 'h::Float' : height of right triangle with hypontenuse, r1, and base, a, colinear with panel.
 - 'a::Float' : length of base of right triangle with height, h, and hypontenuse, r1.

"""
function get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)
    return 1.0 / (2.0 * pi) * (h * (theta2 - theta1) - dmag + a * ln1 - (a - dmag) * ln2)
end

"""
    get_psitildegamma(psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\widetilde{\\Psi}^\\gamma\$

**Arguments:**
 - 'psibargamma::Float' : value of \$\\overbar{\\Psi}^\\gamma\$
 - 'r1mag::Float' : distance from node1 to evaluation point
 - 'r2mag::Float' : distance from node2 to evaluation point
 - 'theta1::Float' : angle between panel and vector from node1 to evaluation point
 - 'theta2::Float' : angle between panel and vector from node2 to evaluation point
 - 'ln1::Float' : value of ln(rmag1), which may be that or 0.0, depening on evaluation point location
 - 'ln2::Float' : value of ln(rmag2), which may be that or 0.0, depening on evaluation point location
 - 'dmag::Float' : panel length
 - 'h::Float' : height of right triangle with hypontenuse, r1, and base, a, colinear with panel.
 - 'a::Float' : length of base of right triangle with height, h, and hypontenuse, r1.
Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
"""
function get_psitildegamma(psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a)
    return a / dmag * psibargamma +
           1 / (4 * pi * dmag) * (r2mag^2 * ln2 - r1mag^2 * ln1 - r2mag^2 / 2 + r1mag^2 / 2)
end

"""
    get_vortex_influence(node1, node2, point)

Calculate vortex influence coefficients on the evaluation point from the panel between node1 and node2.

**Arguments:**
 - 'node1::Array{Float}(2)' : [x y] location of node1
 - 'node2::Array{Float}(2)' : [x y] location of node2
 - 'point::Array{Float}(2)' : [x y] location of evaluation point

"""
function get_vortex_influence(node1, node2, point)

    # Use inputs to get raw distances
    r1, r1mag = get_r(node1, point)
    r2, r2mag = get_r(node2, point)
    d, dmag = get_d(node1, node2)

    #check if point resides on either node and create convenience flag
    if node1 == point
        p1 = true
        p2 = false
    elseif node2 == point
        p1 = false
        p2 = true
    else
        p1 = false
        p2 = false
    end

    # Calculate secondary distances and angles (taking into account whether or not the point lies on one of the nodes)
    if p1 || p2
        a = dmag
        h = 0.0
        if p1
            ln1 = 0.0
            ln2 = log(r2mag)
            theta1 = pi / 2.0 #note that angles don't matter in these cases, since h=0
            theta2 = pi
        else
            ln1 = log(r1mag)
            ln2 = 0.0
            theta1 = 0.0 #note that angles don't matter in these cases, since h=0
            theta2 = pi / 2.0
        end
    else
        theta1 = get_theta(r1, r1mag, d, dmag)
        theta2 = get_theta(r2, r2mag, d, dmag)
        h = get_h(dmag, r1mag, r2mag)
        a = get_a(r1mag, dmag, theta1)
        ln1 = log(r1mag)
        ln2 = log(r2mag)
    end

    # get psibargamma value
    psibargamma = get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)

    # get psitildegamma value
    psitildegamma = get_psitildegamma(
        psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a
    )

    # put psi's together
    return psibargamma - psitildegamma, psitildegamma
end

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
    amat .= 0.0

    # loop through setting up influence coefficients
    for i in 1:N
        for j in 1:(N - 1)

            # obtain influence coefficient for ith evaluation point and j and j+1 panel
            aij, aijp1 = FLOWFoil.get_vortex_influence(nodes[j], nodes[j + 1], nodes[i])

            # add coefficients to matrix at correct nodes
            amat[i, j] += aij
            amat[i, j + 1] += aijp1
        end
    end

    # add in trailing edge contributions
    if meshsystem.meshes[1].blunt_te
        #TODO: Fix this.  Implementation is wrong for Blunt TE
        # get bisection vector
        # get vector along first panel
        d1, _ = get_d(nodes[2], nodes[1])

        # get vector along second panel
        dn, _ = get_d(nodes[end - 1], nodes[end])

        # calculate vector that bisects the first and last panel vectors
        bisector = d1 * LinearAlgebra.norm(dn) + dn * LinearAlgebra.norm(d1)

        # normalize to get the unit vector
        ttehat = bisector / LinearAlgebra.norm(bisector)

        # get panel vector
        dte, _ = get_d(nodes[end], nodes[1])

        # normalize panelvector
        ptehat = dte / LinearAlgebra.norm(dte)

        # get dot product of bisection vector and panel vector.
        tdp = abs(LinearAlgebra.dot(ttehat, ptehat))

        # get cross product of bisection vector and panel vector
        txp = LinearAlgebra.norm(
            LinearAlgebra.cross([ttehat[1]; ttehat[2]; 0.0], [ptehat[1]; ptehat[2]; 0.0])
        )

        # Apply trailing edge contributions
        amat[:, 1] .-= 0.5 * (tdp + txp)
        amat[:, end] .+= 0.5 * (tdp + txp)
    else
        # Replace Nth row of the matrix with the extrapolation of the mean vortex strength to the trailing edge.
        amat[end, :] .= 0.0
        amat[end, 1] = 1.0
        amat[end, 2] = 2.0
        amat[end, 3] = -1.0
        amat[end, end - 2] = 1.0
        amat[end, end - 1] = -2.0
        amat[end, end] = -1.0
    end

    return amat
end

"""
    assemblematrixa!(amat)

Update vortex coefficient matrix to be the full N+1 x N+1 system, including Kutta condition.

**Arguments:**
 - 'amat::Array{Float,2}' : NxN vortex coefficient matrix

"""
function assemblematrixa!(amat)

    # get size of amat:
    r, c = size(amat)

    # Add column of ones to end of matrix for Psi_0
    psi0 = -1.0 * ones(r)
    amat = [amat psi0]

    # add kutta condition row to bottom of matrix
    kutta = zeros(c + 1)
    kutta[1] = 1.0
    kutta[end - 1] = 1.0
    amat = [amat; kutta']

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
function assembleboundaryconditions(meshsystem, freestream)

    # get node locations for convenience
    nodes = meshsystem.meshes[1].airfoil_nodes
    N = length(nodes)

    # Get angle of attack
    alpha = freestream.angleofattack[1] * pi / 180.0

    # Get freestream magnitude
    re = freestream.reynolds[1]
    rho = freestream.density[1]
    mu = freestream.dynamicviscosity[1]
    chord = maximum(getindex.(nodes, 1)) - minimum(getindex.(nodes, 1))
    Vinf = re * mu / (rho * chord)

    # generate boundary condition array
    psi_inf = Vinf * [nodes[i][2] * cos(alpha) - nodes[i][1] * sin(alpha) for i in 1:N]

    # if closed trailing edge, set last element of psi_inf to zero
    if !meshsystem.meshes[1].blunt_te
        psi_inf[end] = 0.0
    end
    # For Kutta Condition, RHS = 0.0
    push!(psi_inf, 0.0)

    return psi_inf
end
