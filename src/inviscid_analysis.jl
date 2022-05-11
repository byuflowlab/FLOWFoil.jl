#=
Inviscid Panel Method Solver

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    function_name(args; kwargs)

Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
"""
function get_psibargamma(theta1, theta2, ln1, ln2, dmag, h, a)
    return 1 / (2 * pi) * (h * (theta2 - theta1) - dmag + a * ln1 - (a - dmag) * ln2)
end

"""
    function_name(args; kwargs)

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
    function_name(args; kwargs)

Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
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
        theta1 = 0.0
        theta2 = 0.0
        if p1
            ln1 = 0.0
            ln2 = log(r2mag)
        else
            ln1 = log(r1mag)
            ln2 = 0.0
        end
    else
        theta1 = get_theta(r1, d)
        theta2 = get_theta(r2, d)
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

"""
    function_name(args; kwargs)

Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
"""
function assemblevortexcoefficients(meshsystem)

    # get system size
    N, Ns = FLOWFoil.sizesystem(meshsystem)

    # get nodes for convenience
    nodes = meshsystem.meshes[1].airfoil_nodes

    #initialize matrix
    a = Array{Float64,2}(undef, N, N)

    # loop through setting up influence coefficients
    for i in 1:N
        for j in 1:N-1

            # obtain influence coefficient for ith evaluation point and j and j+1 panel
            aij, aijp1 = FLOWFoil.get_vortex_influence(nodes[j], nodes[j + 1], nodes[i])

            # add to matrix
            a[i, j] += aij
            a[i, j + 1] += aijp1
        end
    end

    return a
end

"""
    function_name(args; kwargs)

Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
"""
function assemblematrixa!(a)

    # get size of a:
    r, c = size(a)

    # Add column of ones to end of matrix for Psi_0
    psi0 = ones(r)
    a = [a psi0]

    # add kutta condition row to bottom of matrix
    kutta = zeros(c + 1)
    kutta[1] = 1.0
    kutta[end - 1] = 1.0
    a = [a; kutta']

    return a
end

"""
"""
function assemblematrixa(meshsystem)

    #get NxN square of a matrix
    a = assemblevortexcoefficients(meshsystem)

    # get size of a:
    r, c = size(a)

    # Add column of ones to end of matrix for Psi_0
    psi0 = ones(r)
    a = [a psi0]

    # add kutta condition row to bottom of matrix
    kutta = zeros(c + 1)
    kutta[1] = 1.0
    kutta[end - 1] = 1.0
    a = [a; kutta']

    return a
end
"""
    function_name(args; kwargs)

Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
"""
function assembleboundaryconditions(meshsystem, freestream)

    # get node locations for convenience
    nodes = meshsystem.meshes[1].airfoil_nodes
    N = length(nodes)

    # Get angle of attack
    alpha = freestream.anglesofattack[1] * pi / 180.0

    # Get freestream magnitude
    re = freestream.reynolds[1]
    rho = freestream.density[1]
    mu = freestream.dynamicviscosity[1]
    chord = maximum(getindex.(nodes, 1)) - minimum(getindex.(nodes, 1))
    Vinf = re * mu / (rho * chord)

    # generate boundary condition array
    psi_inf = Vinf * [nodes[i][2] * cos(alpha) - nodes[i][1] * sin(alpha) for i in 1:N]
    push!(psi_inf, 0.0)

    return psi_inf
end
