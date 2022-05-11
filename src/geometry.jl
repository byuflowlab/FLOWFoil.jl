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
        numwake = ceil(Int, numnodes / 10 + 10 * wakelength)
    end

    # Get node locations from x,y coordinates
    airfoil_nodes = [[x[i] y[i]] for i in 1:numnodes]

    # initialize the wake panel start point and direction
    wakestart = airfoil_nodes[1]
    wakedir = [1.0; 0.0]

    # initialize blunt_te=false
    blunt_te = false

    # check if open trailing edge
    if airfoil_nodes[1][1] != airfoil_nodes[end][1] ||
        airfoil_nodes[1][2] != airfoil_nodes[end][2]

        # add panel across the trailing edge gap
        push!(airfoil_nodes, airfoil_nodes[1])

        # and update the wake panel starting location to be the midpoint of the gap panel.
        wakestart = (airfoil_nodes[end] .+ airfoil_nodes[end - 1]) / 2.0

        # and update the wake panel initial direction to be the normal of that panel
        wakedir = FLOWFoil.get_normal(airfoil_nodes[end - 1], airfoil_nodes[end])

        # and set blunt_te to true
        blunt_te = true

    else #(closed trailing edge)

        # update wake panel direction to be bisection of trailing edge panel vectors
        # get vector along first panel
        a1 = airfoil_nodes[1] - airfoil_nodes[2]

        # get vector along second panel
        an = airfoil_nodes[end] - airfoil_nodes[end - 1]

        # calculate vector that bisects the first and last panel vectors
        bisector = a1 * LinearAlgebra.norm(an) + an * LinearAlgebra.norm(a1)

        # normalize to get the unit vector
        wakedir = bisector / LinearAlgebra.norm(bisector)
    end

    #TODO: There is something in the method about a trailing half panel, find out what that means and if you should remove the final wake node or not.
    # get initial wake geometry: equidistant panels starting at wakestart point and extending the input percentage of the airfoil chord in the calculated direction.
    wake_nodes = [
        wakestart .+ x .* wakedir for
        x in range(0.0; stop=wakelength * chordlength, length=numwake)
    ]

    # get wake midpoints as well
    wake_midpoints = [
        [(wake_nodes[i + 1][1] + wake_nodes[i][1]) / 2.0 (
            wake_nodes[i + 1][2] + wake_nodes[i][2]
        ) / 2.0] for i in 1:(numwake - 1)
    ]

    # generate mesh object
    mesh = Mesh(airfoil_nodes, wake_nodes, wake_midpoints, blunt_te)

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
    rmag = sqrt((point[1] - node[1])^2 + (point[2] - node[2])^2)

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
    get_theta(r, d)

Get angle (in radians) between panel and vector from node to evaluation point.

**Arguments:**
 - 'r::Vector{Float}' : vector from node to evaluation point
 - 'rmag::Vector{Float}' : distance from node to evaluation point
 - 'd::Vector{Float}' : vector describing panel (vector between adjacent nodes)
 - 'dmag::Vector{Float}' : panel length

"""
function get_theta(r, rmag, d, dmag)

    # use formula for angle between vectors
    # dot product of vectors
    num = LinearAlgebra.dot(r, d)

    # product of magnitude of vectors
    den = rmag * dmag

    # inverse cosine of quotient is angle
    theta = acos(num / den)

    return theta
end

"""
    get_h(dmag, rmag1, rmag2)

Calculate height, h, of triangle formed by panel and vectors from panel nodes to evalulation point.

**Arguments:**
 - 'dmag::Float' : panel length
 - 'rmag1::Float' : distance from node1 to evalulation point
 - 'rmag2::Float' : distance from node2 to evalulation point

"""
function get_h(dmag, rmag1, rmag2)

    # get half perimeter of triangle made between panel and evaluation point
    s = (dmag + rmag1 + rmag2) / 2.0

    # use Heron's formula to find area
    area = sqrt(s * (s - dmag) * (s - rmag1) * (s - rmag2))

    # get h from triangle area formula, where h is the triangle height
    h = 2.0 * area / dmag

    return h
end

"""
    get_a(rmag1, dmag, theta1)

Calculate length of side of triagle colinear with the panel that forms a right triangle with vertex at node 1 and evaluation point.

**Arguments:**
 - 'rmag1::Float' : distance from node1 to evaluation point (hypotenuse)
 - 'dmag::Float' : panel length (colinear with distance, a)
 - 'theta1::Float' : angle between panel and vector from node1 to evaluation point (in radians)

"""
function get_a(rmag1, dmag, theta1)

    # 3 cases depending on size of theta1
    if theta1 > pi / 2.0
        # if angle is greater than 90 degrees, get adjacent angle to small triangle and add to panel length
        a = rmag1 * cos(pi - theta1) + dmag
    elseif theta1 == pi / 2.0
        # if angle is 90 degees, length IS panel length
        a = dmag
    else
        # otherwise, simple trig gives length
        a = rmag1 * cos(theta1)
    end

    return a
end

"""
    get_normal(node1, node2; normalout)

Calculate normal (in or out depending on keyword argument flag) of panel between node1 and node2.

**Arguments:**
 - 'node1::Array{Float}(2)' : [x y] location of node1
 - 'node2::Array{Float}(2)' : [x y] location of node2

**Keyword Arguments:**
 - 'normalout::Bool' : flag whether normal should be out of the body (default = true).

"""
function get_normal(node1, node2; normalout=true)

    # choose 2D rotation matrix depending on whether normal should be in or out of body
    if normalout
        rot = [0.0 -1.0; 1.0 0.0]
    else
        rot = [0.0 1.0; -1.0 0.0]
    end

    # calculate raw normal vector
    normal = rot * (node2' .- node1')
    # get norm of raw vector
    normalnorm = sqrt(normal[1]^2 + normal[2]^2)

    # calculate unit normal vector
    nhat = normal / normalnorm

    # return in correct format
    return [nhat[1] nhat[2]]
end
