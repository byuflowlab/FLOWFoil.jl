#=
Geometry Engine for FLOWFoil.jl

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    generatemesh(x, y; chordlength, wakelength)

Create panels from input geometry coordinates.

**Arguments:**
 - 'x::Vector{Float}' : x coordinates defining airfoil geometry.
 - 'y::Vector{Float}' : y coordinates defining airfoil geometry.

**Keyword Arguments:**
 - chordlength::Float' : length of chord (default = 1)
 - wakelength::Float' : length of wake relative to chord (default = 1)

**Returns**
 - mesh::Mesh : Geometry mesh, including panel nodes, wake nodes, and trailing edge condition.
"""
function generatemesh(x, y; chordlength=1.0, wakelength=1.0)

    # check x and y are equal lengths
    if length(x) != length(y)
        @error("x and y vectors must be of the same length")
    else
        # get number of airfoil nodes for convenience
        numnodes = length(x)
        # get number of wake nodes for convenience
        #  numwake = ceil(Int, numnodes / 10 + 10 * wakelength)
    end

    # Get node locations from x,y coordinates
    airfoil_nodes = [[x[i] y[i]] for i in 1:numnodes]

    # initialize the wake panel start point and direction
    # wakestart = airfoil_nodes[1]
    # wakedir = [1.0; 0.0]

    # initialize blunt_te=false
    blunt_te = false

    # check if open trailing edge
    if airfoil_nodes[1][1] != airfoil_nodes[end][1] ||
        airfoil_nodes[1][2] != airfoil_nodes[end][2]

        #TODO: probably not...
        # add panel across the trailing edge gap
        #        push!(airfoil_nodes, airfoil_nodes[1])

        #TODO LATER (probably in different function)
        # update the wake panel starting location to be the midpoint of the gap panel.
        #        wakestart = (airfoil_nodes[end] .+ airfoil_nodes[end - 1]) / 2.0

        # and update the wake panel initial direction to be the normal of that panel
        #       wakedir = FLOWFoil.get_normal(airfoil_nodes[end - 1], airfoil_nodes[end])

        # set blunt_te to true
        blunt_te = true

    else #(closed trailing edge)

        #TODO: probably put this elsewhere
        # update wake panel direction to be bisection of trailing edge panel vectors
        # get vector along first panel
        #        a1 = airfoil_nodes[1] - airfoil_nodes[2]

        # get vector along second panel
        #       an = airfoil_nodes[end] - airfoil_nodes[end - 1]

        # calculate vector that bisects the first and last panel vectors
        #      bisector = a1 * LinearAlgebra.norm(an) + an * LinearAlgebra.norm(a1)

        # normalize to get the unit vector
        #     wakedir = bisector / LinearAlgebra.norm(bisector)
    end

    #TODO: There is something in the method about a trailing half panel, find out what that means and if you should remove the final wake node or not.
    # get initial wake geometry: equidistant panels starting at wakestart point and extending the input percentage of the airfoil chord in the calculated direction.
    # wake_nodes = [
    #    wakestart .+ x .* wakedir for
    #     x in range(0.0; stop=wakelength * chordlength, length=numwake)
    # ]

    # get wake midpoints as well
    #wake_midpoints = [
    #    [(wake_nodes[i + 1][1] + wake_nodes[i][1]) / 2.0 (
    #        wake_nodes[i + 1][2] + wake_nodes[i][2]
    #    ) / 2.0] for i in 1:(numwake - 1)
    #]

    # generate mesh object
    mesh = FLOWFoil.Mesh(airfoil_nodes, blunt_te)

    return mesh
end

"""
    sizesystem(meshsystem)

Count size of inviscid system matrix.

**Arguments:**
 - 'meshsystem::MeshSystem' : The system for which to calculate the linear system size.
"""
function sizesystem(meshsystem)

    # initialize
    # number of bodies for convenience
    numbodies = length(meshsystem.meshes)

    # initialize total system size
    N = 0

    # initialize system size contributions from each mesh
    Ns = ones(Int, numbodies)

    # Count number of airfoil nodes in each mesh.
    for i in 1:numbodies
        Ns[i] = length(meshsystem.meshes[i].airfoil_nodes)
        N += Ns[i]
    end

    return N, Ns
end

"""
    function get_r(node,point)

Calculate the vector, \$\\mathbf{r}\$, and distance, \$|r|\$, from the node to the evaluation point

**Arguments:**
 - `node::Array{Float}` : [x y] position of node
 - `point::Array{Float}` : [x y] position of point.

**Returns**
 - 'r::Vector{Float}' : vector from node to evaluation point
 - 'rmag::Float' : length of panel between node and evaluation point
"""
function get_r(node, point)

    # Calculate vector
    r = point .- node

    # Calculate magnitude
    rmag = sqrt(r[1]^2 + r[2]^2)

    return r, rmag
end

"""
    get_d(node1, node2)

Calculate panel length (between adjacent nodes).


**Arguments:**
 - 'node1::Array{Float}(2)' : [x y] location of first node
 - 'node2::Array{Float}(2)' : [x y] location of second node

**Returns**
 - 'd::Vector{Float}' : vector from node1 to node2
 - 'dmag::Float' : length of panel between node1 and node2
"""
function get_d(node1, node2)

    # simply call get_r, since it's exactly what is needed
    return get_r(node1, node2)
end

"""
    get_theta(h, a)

Get angle (in radians) between panel and vector from node to evaluation point.

**Arguments:**
 - 'h::Float' : Distance, normal to panel, between panel and evaluation point.
 - 'a::Float' : Distance, tangent to panel, between node1 and evaluation point.

"""
function get_theta(h, a)
    return atan(h, a)
end

"""
    get_theta(h, a, dmag)

Get angle (in radians) between panel and vector from node to evaluation point.

**Arguments:**
 - 'h::Float' : Distance, normal to panel, between panel and evaluation point.
 - 'a::Float' : Distance, tangent to panel, between node1 and evaluation point.
 - 'dmag::Float' : Panel lentgh.

"""
function get_theta(h, a, dmag)
    return atan(h, a - dmag)
end

"""
    get_h(r1, d, dmag)

Calculate distance from panel to evalulation point in the panel normal direction.

**Arguments:**
 - 'r1::Vector{Float}' : vector from node1 to evalulation point.
 - 'd::Vector{Float}' : vector from node1 to node2.
 - 'dmag::Float' : panel length

"""
function get_h(r1, d, dmag)

    # get unit normal to panel
    nhat = get_normal(d, dmag)

    # calculate h (dot product of unit normal and r1 vector
    h = r1[1] * nhat[1] + r1[2] * nhat[2]

    return h
end

"""
    get_a(r1, d, dmag)

Calculate distance from panel to evalulation point in the panel tangent direction.

**Arguments:**
 - 'r1::Vector{Float}' : vector from node1 to evalulation point.
 - 'd::Vector{Float}' : vector from node1 to node2.
 - 'dmag::Float' : panel length
    get_a(rmag1, dmag, theta1)

"""
function get_a(r1, d, dmag)

    # Get unit tangent vector
    that = get_tangent(d, dmag)

    # calculate a (dot product of unit tangent and r1 vector)
    a = r1[1] * that[1] + r1[2] * that[2]

    return a
end

"""
    get_tangent(d, dmag)

Get unit tangent to panel.

**Arguments:**
 - 'd::Vector{Float}' : vector from node1 to node2.
 - 'dmag::Float' : panel length

"""
function get_tangent(d, dmag)
    return d / dmag
end

"""
    get_normal(d, dmag)

Get unit normal to panel.

**Arguments:**
 - 'd::Vector{Float}' : vector from node1 to node2.
 - 'dmag::Float' : panel length

"""
function get_normal(d, dmag)

    # get unit tangent
    that = get_tangent(d, dmag)

    # use fancy trick to rotate to be unit normal
    nhat = [-that[2]; that[1]]

    return nhat
end

"""
    get_distances(node1, node2, point)

Get vectors and magnitudes for panel and between nodes and validation points.

**Arguments:**
 - 'node1::Array{Float}' : [x y] position of node1.
 - 'node2::Array{Float}' : [x y] position of node2.
 - 'point::Array{Float}' : [x y] position of evaluation point.

**Returns:**
 - 'r1::Vector{Float}' : vector from node1 to evaluation point.
 - 'r1mag::Float' : distance from node1 to evaluation point.
 - 'r2::Vector{Float}' : vector from node2 to evaluation point.
 - 'r2mag::Float' : distance from node2 to evaluation point.
 - 'd::Vector{Float}' : vector from node1 to node2.
 - 'dmag::Float' : panel length

"""
function get_distances(node1, node2, point)
    r1, r1mag = get_r(node1, point)
    r2, r2mag = get_r(node2, point)
    d, dmag = get_d(node1, node2)
    return r1, r1mag, r2, r2mag, d, dmag
end

"""
    get_orientation(node1, node2, point)

Get angles between panel and evaluation point, ln of distances from nodes to evaluation point, and evaluation point position relative to panel.

**Arguments:**
 - 'node1::Array{Float}' : [x y] position of node1.
 - 'node2::Array{Float}' : [x y] position of node2.
 - 'point::Array{Float}' : [x y] position of evaluation point.

**Returns:**
 - 'theta1::Float' : Angle between panel and evaluation point, centered at node1.
 - 'theta2::Float' : Angle between panel and evaluation point, centered at node2.
 - 'ln1::Float' : Natural log of distance from node1 to evaluation point.
 - 'ln2::Float' : Natural log of distance from node2 to evaluation point.
 - 'h::Float' : Distance from panel to evaluation in panel normal direction.
 - 'a::Float' : Distance from node1 to evaluation in panel tangent direction.

"""
function get_orientation(node1, node2, point)

    # get distances
    r1, r1mag, r2, r2mag, d, dmag = get_distances(node1, node2, point)

    # Get distances normal and tangent to panel from node1
    h = get_h(r1, d, dmag)
    a = get_a(r1, d, dmag)

    #check if point resides on either node and create convenience flag
    #TODO: probably want to set this up with haveing some sort of tolerance rather than an exact value in case user data is not exact, but close.
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
    if p1
        ln1 = 0.0
        ln2 = log(r2mag)
        theta1 = pi
        theta2 = pi
    elseif p2
        ln1 = log(r1mag)
        ln2 = 0.0
        theta1 = 0.0
        theta2 = 0.0
    else
        ln1 = log(r1mag)
        ln2 = log(r2mag)
        theta1 = get_theta(h, a)
        theta2 = get_theta(h, a, dmag)
    end

    return theta1, theta2, ln1, ln2, h, a
end
