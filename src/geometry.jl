#=
Geometry Engine for FLOWFoil.jl

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
=#

"""
    generatemesh(coordinates; order="linear")

Create panels from input geometry coordinates.

**Arguments:**
 - 'x::Vector{Float}' : x coordinates defining airfoil geometry.
 - 'y::Vector{Float}' : y coordinates defining airfoil geometry.

**Returns**
 - mesh::Mesh : Geometry mesh, including panel edge points and collocation points.
"""
function generatemesh(x, y; chordlength=1.0, wakelength=1.0)

    # check x and y are equal lengths
    if length(x) != length(y)
        @error("x and y vectors must be of the same length")
    else
        # get number of nodes for convenience
        numnodes = length(x)
        # get number of wake nodes for convenience
        numwake = ceil(Int, numnodes / 10 + 10 * wakelength)
    end

    # Get panel edges from x,y coordinates
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
        wakestart = (airfoil_nodes[end] .+ airfoil_nodes[end - 1])/2.0

        # and update the wake panel initial direction to be the normal of that panel
        wakedir = FLOWFoil.get_normal(airfoil_nodes[end - 1], airfoil_nodes[end])

        # and set blunt_te to true
        blunt_te = true

    else

        # update wake panel direction to be mean of trailing edge panel angles

        # get vector along first panel
        a1 = airfoil_nodes[1] - airfoil_nodes[2]

        # get vector along second panel
        an = airfoil_nodes[end] - airfoil_nodes[end - 1]

        # calculate vector that bisects the first and last panel vectors
        bisector = a1 * LinearAlgebra.norm(an) + an * LinearAlgebra.norm(a1)

        # normalize to get the unit vector
        wakedir = bisector / LinearAlgebra.norm(bisector)
    end

    # get initial wake geometry: equidistant panels starting at wakestart point and extending the input percentage of the airfoil chord in the calculated direction.
    wake_nodes = [        wakestart .+ x .* wakedir  for        x in range(0.0; stop=wakelength * chordlength, length=numwake)    ]

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
    # number of panels for convenience
    numpanels = length(meshsystem.meshes)

    # total system size
    N = 0

    # system size contributions from each mesh
    Ns = ones(Int, numpanels)

    # Count number of collocation points in each mesh.
    for i in 1:numpanels
        Ns[i] = length(meshsystem.meshes[i].airfoil_nodes)
        N += Ns[i]
    end

    return N, Ns
end

"""
    function get_r(edge,point)

Calculate the vector, \$\\mathbf{r}\$, and distance, \$|r|\$, from the node to the evaluation point

**Arguments:**
 - `node::Array{Float}` : [x y] position of node
 - `point::Array{Float}` : [x y] position of point.
"""
function get_r(node, point)

    # Calculate vector
    r = point .- node

    # Calculate magnitude
    magr = sqrt((point[1] - node[1])^2 + (point[2] - node[2])^2)

    return r, magr
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
function get_d(node1, node2)

    # simply call get_r, since it's exactly what is needed
    return get_r(node1, node2)
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
function get_theta(r, d)

    # use formula for angle between vectors
    # dot product of vectors
    num = LinearAlgebra.dot(r, d)

    # product of magnitude of vectors
    den = LinearAlgebra.norm(r) * LinearAlgebra.norm(d)

    # inverse cosine of quotient is angle
    theta = acos(num / den)

    return theta
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
    function_name(args; kwargs)

Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
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
    function_name(args; kwargs)

Function Description.

Detailed Description.

**Arguments:**
 - arg::type : description.

**Returns**
 - output::type : description.
"""
function get_normal(e1, e2; normalout=true)
    if normalout
        rot = [0.0 -1.0; 1.0 0.0]
    else
        rot = [0.0 1.0; -1.0 0.0]
    end

    nhat = rot * (e2' .- e1')
    nhatnorm = sqrt(nhat[1]^2 + nhat[2]^2)

    normal = nhat / nhatnorm

    return [normal[1] normal[2]]
end

"""
    distances(meshsystem)
 - meshsystem::MeshSystem : Mesh System for which to calculate distances.

**Returns**
 - rs::Array{Float,2} : Distances from panel edges to control points.
"""
function distances(meshsystem)

    # get reqired matrix size
    N, Ns = FLOWFoil.sizesystem(meshsystem)

    # Get number of meshes for convenience
    meshes = meshsystem.meshes
    nm = length(meshes)

    # Initialize output
    # note that there are n+1 edges for each mesh (where n=number of collocation points), so there are N+length(Ns) rows in the output matrices.
    rs = Array{Array{Float64}}(undef, N, N + length(Ns))
    rmags = Array{Float64,2}(undef, N, N + length(Ns))

    # Calculate Distances
    # Loop over meshes (edges)
    for m in 1:nm

        #rename for convenience
        edges = meshes[m].edges

        #get offset location for global matrix: total number of edge points from the meshes thus far.
        joffset = sum(Ns[1:m]) - Ns[m] + m - 1

        # Loop over meshes (collocation points)
        for n in 1:nm

            #rename for convenience
            points = meshes[n].collocation_points

            #get offset location for global matrix: total number of collocation points from the meshes thus far.
            ioffset = sum(Ns[1:n]) - Ns[n]

            for j in 1:(Ns[m] + 1)
                for i in 1:Ns[n]
                    # get vector and magnitude
                    vec, mag = FLOWFoil.get_r(edges[j], points[i])

                    # assign vector to global vector matrix
                    rs[i + ioffset, j + joffset] = vec

                    # assign magnitude to global magnitude matrix.
                    rmags[i + ioffset, j + joffset] = mag
                end
            end
        end
    end

    return rs, rmags
end

"""
"""
function normals(meshsystem; normalout=true)

    # get reqired matrix size
    N, Ns = FLOWFoil.sizesystem(meshsystem)

    # Get number of meshes for convenience
    meshes = meshsystem.meshes
    nm = length(meshes)

    # Initialize output
    panelnormals = Array{Array{Float64}}(undef, N)

    # Loop over meshes
    for n in 1:nm

        #rename for convenience
        edges = meshes[n].edges

        #get offset location for global matrix: total number of panels from the meshes thus far.
        ioffset = sum(Ns[1:n]) - Ns[n]

        for i in 1:Ns[n]

            # get unit normal
            panelnormals[i + ioffset] = FLOWFoil.get_normal(
                edges[i], edges[i + 1]; normalout=normalout
            )
        end
    end

    return panelnormals
end

