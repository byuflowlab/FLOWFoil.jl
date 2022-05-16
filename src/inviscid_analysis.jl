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
    return 1.0 / (2.0 * pi * dmag) *
           (h * (theta2 - theta1) - dmag + a * ln1 - (a - dmag) * ln2)
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

"""
function get_psitildegamma(psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a)
    return a * psibargamma +
           1 / (4 * pi * dmag) * (r2mag^2 * ln2 - r1mag^2 * ln1 - r2mag^2 / 2 + r1mag^2 / 2)
end

"""
    get_psibarsigma(theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\overbar{\\Psi}^\\sigma\$

**Arguments:**
 - 'theta1::Float' : Angle between panel and evaluation point, centered at node1.
 - 'theta2::Float' : Angle between panel and evaluation point, centered at node2.
 - 'ln1::Float' : Natural log of distance from node1 to evaluation point.
 - 'ln2::Float' : Natural log of distance from node2 to evaluation point.
 - 'h::Float' : Distance from panel to evaluation in panel normal direction.
 - 'a::Float' : Distance from node1 to evaluation in panel tangent direction.
"""
function get_psibarsigma(theta1, theta2, ln1, ln2, dmag, h, a)
    return 1 / (2 * pi) * (a * (theta1 - theta2) + dmag * theta2 + h * ln1 - h * ln2)
end

"""
    get_psitildesigma(psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a)

Calculate value of  \$\\widetilde{\\Psi}^\\sigma\$

**Arguments:**
 - 'psibargamma::Float' : value of \$\\overbar{\\Psi}^\\gamma\$
 - 'r1mag::Float' : distance from node1 to evaluation point
 - 'r2mag::Float' : distance from node2 to evaluation point
 - 'theta1::Float' : angle between panel and vector from node1 to evaluation point
 - 'theta2::Float' : angle between panel and vector from node2 to evaluation point
 - 'dmag::Float' : panel length
 - 'h::Float' : height of right triangle with hypontenuse, r1, and base, a, colinear with panel.
 - 'a::Float' : length of base of right triangle with height, h, and hypontenuse, r1.
"""
function get_psitildesigma(psibarsigma, r1mag, r2mag, theta1, theta2, dmag, h, a)
    return a / dmag * phibarsigma +
           1.0 / (4 * pi * dmag) * (rmag2^2 * theta2 - rmag1^2 * theta1 - h * dmag)
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
    r1, r1mag, r2, r2mag, d, dmag = get_distances(node1, node2, point)

    # Calculate a, h, and natural logs based on position of point
    theta1, theta2, ln1, ln2, h, a = get_orientation(node1, node2, point)

    # get psibargamma value
    psibargamma = get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)

    # get psitildegamma value
    psitildegamma = get_psitildegamma(
        psibargamma, r1mag, r2mag, theta1, theta2, ln1, ln2, dmag, h, a
    )

    # put psi's together
    return (psibargamma - psitildegamma), psitildegamma
end

"""
    get_source_influence(node1, node2, point)

Calculate source influence coefficients on the evaluation point from the panel between node1 and node2.

**Arguments:**
 - 'node1::Array{Float}(2)' : [x y] location of node1
 - 'node2::Array{Float}(2)' : [x y] location of node2
 - 'point::Array{Float}(2)' : [x y] location of evaluation point
"""
function get_source_influence(node1, node2, point)

    # Use inputs to get raw distances
    r1, r1mag, r2, r2mag, d, dmag = get_distances(node1, node2, point)

    # Calculate a, h, and natural logs based on position of point
    theta1, theta2, ln1, ln2, h, a = get_orientation(node1, node2, point)

    #get psibarsigma value
    psibarsigma = get_psibarsigma(theta1, theta2, ln1, ln2, dmag, h, a)

    # shift source in order to get a better behaved branch cut orientation
    if (theta1 + theta2) > pi
        psibarsigma -= 0.25 * dmag
    else
        psibarsigma += 0.75 * dmag
    end

    return psibarsigma
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

        # get bisection vector
        # get vector along first panel
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
        amat[end, end] = 1.0
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
    psi_inf = Vinf * [-nodes[i][2] * cos(alpha) + nodes[i][1] * sin(alpha) for i in 1:N]

    # if closed trailing edge, set last element of psi_inf to zero
    if !meshsystem.meshes[1].blunt_te
        psi_inf[end] = 0.0
    end
    # For Kutta Condition, RHS = 0.0
    push!(psi_inf, 0.0)

    return psi_inf
end
