#=
Geometry Engine for FLOWFoil.jl

Authors: Judd Mehr,

Date Started: 27 April 2022

Change Log:
10/22 - Add axisymmetric geometry functions
=#

"""
    generate_mesh(x, y; chordlength, wakelength)

Create panels from input geometry coordinates.

**Arguments:**
 - `x::Vector{Float}` : x coordinates defining airfoil geometry.
 - `y::Vector{Float}` : y coordinates defining airfoil geometry.

**Keyword Arguments:**
 - `gaptolerance::Float` : Tolerance for how close, relative to the chord, the trailing edge nodes can be before being considered a sharp trailing edge. (default = 1e-10)

**Returns**
 - `mesh::PlanarMesh` : Geometry mesh, including panel nodes and trailing edge condition.
"""
function generate_mesh(x, y; gaptolerance=1e-10)

    # check x and y are equal lengths
    if length(x) != length(y)
        @error("x and y vectors must be of the same length")
    else
        # get number of airfoil nodes for convenience
        numnodes = length(x)
    end

    # Get node locations from x,y coordinates
    nodes = [[x[i] y[i]] for i in 1:numnodes]

    # get trailing edge information
    tdp, txp, trailing_edge_gap = get_trailing_edge_info(nodes)

    # get chord length
    chordlength = maximum(x) - minimum(x)

    # check if open trailing edge
    if abs(trailing_edge_gap) > gaptolerance * chordlength

        # set blunt_te to true
        blunt_te = true

    else #(closed trailing edge)

        #set blunt_te to false
        blunt_te = false
    end

    # generate mesh object
    mesh = FLOWFoil.PlanarMesh(nodes, chordlength, blunt_te, trailing_edge_gap, tdp, txp)

    return mesh
end

"""
    generate_mesh(coordinates; kwargs)

Identical to implementation with x and y separate, but here with x,y coordinates together in a single array [X Y].

**Arguments:**
 - `coordinates::Array{Float,2}` : array of both x and y coordinates (x first column, y second column).
"""
function generate_mesh(coordinates; gaptolerance=1e-10)

    # Separate out coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    return generate_mesh(x, y; gaptolerance=gaptolerance)
end

"""
    position_coordinates!(meshes, scales, angles, locations)

Take in meshes and adjust scale, leading edge location, and angle of attack of the individual meshes in the system.  Updates mesh objects in place.

**Arguments:**
 - `meshes::Array{PlanarMesh}` : Array of mesh objects.
 - `scales::Array{Float}` : Array of numbers by which to scale respective meshes.
 - `angles::Array{Float}` : Array of angles, in degrees, by which to rotate respective meshes (positive = pitch up).
 - `locations::Array{Array{Float}}` : Array of [x y] positions of leading edges for respective meshes.

"""
function position_coordinates(coordinates, scale, angle, location; flipped=false)

    #flip if needed
    if flipped
        coordinates[:, 2] .*= -1.0
        reverse!(coordinates; dims=1)
    end

    # scale
    coordinates .*= scale

    # get rotation matrix
    R = [cosd(-angle) -sind(-angle); sind(-angle) cosd(-angle)]

    # rotate and translate
    for j in 1:length(coordinates[:, 1])
        coordinates[j, :] = R * coordinates[j, :]
        coordinates[j, :] .+= location
    end

    return coordinates[:, 1], coordinates[:, 2]
end

"""
    position_coordinates!(meshsystem)

Identical to position_coordinates!, but taking the inputs in as a PlanarMeshSystem object.

**Arguments:**
 - `meshsystem::PlanarMeshSystem` : Mesh system object to position.
"""
function position_coordinates!(meshsystem)
    position_coordinates!(
        meshsystem.meshes, meshsystem.scales, meshsystem.angles, meshsystem.locations
    )

    return nothing
end

"""
    size_system(meshsystem)

Count size of inviscid system matrix.

**Arguments:**
 - `meshsystem::Array{PlanarMesh}` : The system for which to calculate the linear system size.
"""
function size_system(meshes; axisymmetric=false)

    # initialize
    # number of bodies for convenience
    numbodies = length(meshes)

    # initialize total system size
    N = 0

    # initialize system size contributions from each mesh
    Ns = [1 for i in 1:numbodies]

    # Count number of airfoil nodes in each mesh.
    for i in 1:numbodies
        if axisymmetric
            Ns[i] = length(meshes[i].panels)
        else
            Ns[i] = length(meshes[i].nodes)
        end
        N += Ns[i]
    end

    return N, Ns
end

"""
    get_offset(Ns)

Get the offset values for the mesh system to be used in the system matrix assembly.

**Arguments:**
 - `Ns::Array{Float}` : Array of numbers of nodes for each airfoil in the system.
"""
function get_offset(Ns)
    return [0; cumsum(Ns[1:(end - 1)])]
end

"""
    get_trailing_edge_info(nodes)

Calculate various items needed for trailing edge treatment.

**Arguments:**
 - `nodes::Array{Float,2}` : Array of [x y] locations for the airfoil nodes.

**Returns:**
 - `tdp::Float` : dot product of TE bisection and TE gap unit vectors
 - `txp::Float` : "cross product" of TE bisection and TE gap unit vectors
 - `trailing_edge_gap::Float` : TE gap distance
"""
function get_trailing_edge_info(nodes)
    # get bisection vector
    # get vector along first panel
    d1, d1mag = get_d(nodes[2], nodes[1])
    d1 /= d1mag

    # get vector along second panel
    dn, dnmag = get_d(nodes[end - 1], nodes[end])
    dn /= dnmag

    # calculate vector that bisects the first and last panel vectors
    bisector = 0.5 * (d1 + dn)

    # normalize to get the unit vector
    ttehat = bisector / sqrt(bisector[1]^2 + bisector[2]^2)

    # get panel vector
    dte, dtemag = get_d(nodes[1], nodes[end])

    # normalize panelvector
    dtehat = dte / dtemag

    # get dot product of bisection vector and panel vector.
    tdp = ttehat[1] * dtehat[1] + ttehat[2] * dtehat[2]

    # get cross product of bisection vector and panel vector
    txp = abs(ttehat[1] * dtehat[2] - ttehat[2] * dtehat[1])

    # get trailing edge gap
    trailing_edge_gap = -dte[1] * ttehat[2] + dte[2] * ttehat[1]

    return isnan(tdp) ? 0.0 : tdp, isnan(txp) ? 0.0 : txp, trailing_edge_gap
end

"""
    function get_r(node,point)

Calculate the vector, \$\\mathbf{r}\$, and distance, \$|r|\$, from the node to the evaluation point

**Arguments:**
 - `node::Array{Float}` : [x y] position of node
 - `point::Array{Float}` : [x y] position of point.

**Returns**
 - `r::Vector{Float}` : vector from node to evaluation point
 - `rmag::Float` : length of panel between node and evaluation point
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
 - `node1::Array{Float}(2)` : [x y] location of first node
 - `node2::Array{Float}(2)` : [x y] location of second node

**Returns**
 - `d::Vector{Float}` : vector from node1 to node2
 - `dmag::Float` : length of panel between node1 and node2
"""
function get_d(node1, node2)

    # simply call get_r, since it`s exactly what is needed
    return get_r(node1, node2)
end

"""
    get_theta(h, a)

Get angle (in radians) between panel and vector from node1 to evaluation point.

**Arguments:**
 - `h::Float` : Distance, normal to panel, between panel and evaluation point.
 - `a::Float` : Distance, tangent to panel, between node1 and evaluation point.

"""
function get_theta(h, a)
    return atan(h, a)
end

"""
    get_theta(h, a, dmag)

Get angle (in radians) between panel and vector from node2 to evaluation point.

**Arguments:**
 - `h::Float` : Distance, normal to panel, between panel and evaluation point.
 - `a::Float` : Distance, tangent to panel, between node1 and evaluation point.
 - `dmag::Float` : Panel lentgh.

"""
function get_theta(h, a, dmag)
    return atan(h, a - dmag)
end

"""
    get_h(r1, d, dmag)

Calculate distance from panel to evalulation point in the panel normal direction.

**Arguments:**
 - `r1::Vector{Float}` : vector from node1 to evalulation point.
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

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
 - `r1::Vector{Float}` : vector from node1 to evalulation point.
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

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
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_tangent(d, dmag)
    return (dmag == 0.0) ? [0.0; 0.0] : (d / dmag)
end

"""
    get_normal(d, dmag)

Get unit normal to panel.

**Arguments:**
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

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
 - `node1::Array{Float}` : [x y] position of node1.
 - `node2::Array{Float}` : [x y] position of node2.
 - `point::Array{Float}` : [x y] position of evaluation point.

**Returns:**
 - `r1::Vector{Float}` : vector from node1 to evaluation point.
 - `r1mag::Float` : distance from node1 to evaluation point.
 - `r2::Vector{Float}` : vector from node2 to evaluation point.
 - `r2mag::Float` : distance from node2 to evaluation point.
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_distances(node1, node2, point)

    # Get vector and magnitude from node1 to point
    r1, r1mag = get_r(node1, point)

    # Get vector and magnitude from node2 to point
    r2, r2mag = get_r(node2, point)

    # Get vector and magnitude of panel
    d, dmag = get_d(node1, node2)

    return r1, r1mag, r2, r2mag, d, dmag
end

"""
    get_orientation(node1, node2, point)

Get angles between panel and evaluation point, ln of distances from nodes to evaluation point, and evaluation point position relative to panel.

**Arguments:**
 - `node1::Array{Float}` : [x y] position of node1.
 - `node2::Array{Float}` : [x y] position of node2.
 - `point::Array{Float}` : [x y] position of evaluation point.

**Returns:**
 - `theta1::Float` : Angle between panel and evaluation point, centered at node1.
 - `theta2::Float` : Angle between panel and evaluation point, centered at node2.
 - `ln1::Float` : Natural log of distance from node1 to evaluation point.
 - `ln2::Float` : Natural log of distance from node2 to evaluation point.
 - `h::Float` : Distance from panel to evaluation in panel normal direction.
 - `a::Float` : Distance from node1 to evaluation in panel tangent direction.

"""
function get_orientation(node1, node2, point; epsilon=1e-9)

    # get distances
    r1, r1mag, r2, r2mag, d, dmag = get_distances(node1, node2, point)

    # Get distances normal and tangent to panel from node1
    h = get_h(r1, d, dmag)
    a = get_a(r1, d, dmag)

    # Calculate natural log values, setting to zero if the point is close enough to the node, also get theta values and change based on point location if necessary
    theta1 = get_theta(h, a)
    theta2 = get_theta(h, a, dmag)
    if r1mag < epsilon
        ln1 = 0.0
        theta1 = pi
        theta2 = pi
    else
        ln1 = log(r1mag)
    end

    if r2mag < epsilon
        ln2 = 0.0
        theta1 = 0.0
        theta2 = 0.0
    else
        ln2 = log(r2mag)
    end

    return theta1, theta2, ln1, ln2, h, a
end

#"""
#"""
#function initialize_wake()

#    # initialize the wake panel start point and direction
#    # wakestart = nodes[1]
#    # wakedir = [1.0; 0.0]

#    #TODO LATER (probably in different function)
#    # update the wake panel starting location to be the midpoint of the gap panel.
#    #        wakestart = (nodes[end] .+ nodes[end - 1]) / 2.0

#    # and update the wake panel initial direction to be the normal of that panel
#    #       wakedir = FLOWFoil.get_normal(nodes[end - 1], nodes[end])
#    #TODO: probably put this elsewhere
#    # update wake panel direction to be bisection of trailing edge panel vectors
#    # get vector along first panel
#    #        a1 = nodes[1] - nodes[2]

#    # get vector along second panel
#    #       an = nodes[end] - nodes[end - 1]

#    # calculate vector that bisects the first and last panel vectors
#    #      bisector = a1 * LinearAlgebra.norm(an) + an * LinearAlgebra.norm(a1)

#    # normalize to get the unit vector
#    #     wakedir = bisector / LinearAlgebra.norm(bisector)

#    #TODO: There is something in the method about a trailing half panel, find out what that means and if you should remove the final wake node or not.
#    # get initial wake geometry: equidistant panels starting at wakestart point and extending the input percentage of the airfoil chord in the calculated direction.
#    # wake_nodes = [
#    #    wakestart .+ x .* wakedir for
#    #     x in range(0.0; stop=wakelength * chordlength, length=numwake)
#    # ]

#    # get wake midpoints as well
#    #wake_midpoints = [
#    #    [(wake_nodes[i + 1][1] + wake_nodes[i][1]) / 2.0 (
#    #        wake_nodes[i + 1][2] + wake_nodes[i][2]
#    #    ) / 2.0] for i in 1:(numwake - 1)
#    #]

#end

####################################
##### ----- AXISYMMETRIC ----- #####
####################################

"""
    generate_axisym_mesh(x, r; bodyofrevolution)

Generate mesh for axisymmetric body.

**Arguments:**
- `x::Array{Float}` : x-coordinates of geometry
- `r::Array{Float}` : r-coordinates of geometry

**Keyword Arguments:**
- `bodyofrevolution::Bool` : flag whether body is a body of revolution (default=true)

**Returns:**
- `mesh::FLOWFoil.AxiSymMesh` : axisymmetric mesh object
"""
function generate_axisym_mesh(x, r; bodyofrevolution=true, ex=1e-5)

    #check of any r coordinates are negative
    @assert all(x -> x >= -eps(), r)

    #initialize panels
    panels = Array{AxiSymPanel}(undef, length(x) - 1)

    cpx = [0.0 for i in 1:(length(x) - 1)]
    cpr = [0.0 for i in 1:(length(x) - 1)]
    nhat = [[0.0; 0.0] for i in 1:(length(x) - 1)]
    dmag = [0.0 for i in 1:(length(x) - 1)]
    sine = [0.0 for i in 1:(length(x) - 1)]
    cosine = [0.0 for i in 1:(length(x) - 1)]
    slope = [0.0 for i in 1:(length(x) - 1)]
    curve = [0.0 for i in 1:(length(x) - 1)]

    for i in 1:(length(x) - 1)

        #calculate control point
        cpx[i] = 0.5 * (x[i] + x[i + 1])
        cpr[i] = 0.5 * (r[i] + r[i + 1])

        #calculate length
        d, dmag[i] = get_d([x[i]; r[i]], [x[i + 1]; r[i + 1]])

        #calculate normal
        nhat[i] = get_normal(d, dmag[i])

        #find minimum x point (i.e. the leading edge point
        _, minx = findmin(x)

        #use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        beta = atan(d[2] / d[1])

        #apply corrections as needed based on orientation of panel in coordinate frame.
        if (d[1] < 0.0) && (i > minx)
            #if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            slope[i] = beta - pi

        elseif (d[1] < 0.0) && (i < minx)
            #if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            slope[i] = beta + pi
        else
            slope[i] = beta
        end

        #TODO: This is the version from the book code.  Maybe it's more robust?
        ## sine[i] = (r[i + 1] - r[i]) / dmag[i]
        #sine[i] = (d[2]) / dmag[i]
        #cosine[i] = (d[1]) / dmag[i]
        ## cosine[i] = (x[i + 1] - x[i]) / dmag[i]
        #abscos = abs(cosine[i])
        #if abscos > ex
        #    #use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        #    beta = atan(sine[i] / cosine[i])
        #end

        ##if the panel is nearly vertical, set the panel slope to vertical in the correct direction.
        #if abscos < ex
        #    slope[i] = sign(sine[i]) * pi / 2.0
        #end

        ##otherwise (in most cases)
        #if cosine[i] > ex
        #    slope[i] = beta
        #end

        ## For special cases
        ##find minimum x point (i.e. the leading edge point
        #_, minx = findmin(x)

        ##if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
        #if (cosine[i] < -ex) && (i > minx)
        #    slope[i] = beta - pi
        #end

        ##if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
        #if (cosine[i] < -ex) && (i < minx)
        #    slope[i] = beta + pi
        #end

    end

    for i in 2:(length(x) - 2)
        curve[i] = (slope[i + 1] - slope[i - 1]) / 8.0 / pi
    end

    for i in 1:(length(x) - 1)
        #generate panel objects
        panels[i] = AxiSymPanel([cpx[i]; cpr[i]], dmag[i], nhat[i], slope[i], curve[i])
    end

    return AxiSymMesh(panels, bodyofrevolution)
end

"""
    get_ring_geometry(paneli, panelj)

Obtain relevant geometry associated with ring singularity influence calculations.

**Arguments:**
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).

**Returns:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `dmagj::Float` : length of the jth panel
- `m::Float` : Elliptic Function parameter
- `nhati::Array{Float}` : unit normal vector of ith panel
"""
function get_ring_geometry(paneli, panelj)

    #rename for convenience
    dmagj = panelj.length

    nhati = paneli.normal

    xi = paneli.controlpoint[1]
    ri = paneli.controlpoint[2]

    xj = panelj.controlpoint[1]
    rj = panelj.controlpoint[2]

    #get x and r for these panels
    x = (xi - xj) / rj
    r = ri / rj

    #get phi for these panels
    m = 4.0 * r / (x^2 + (r + 1.0)^2)
    if isnan(m)
        println("m: ", m)
        display(paneli.controlpoint)
        display(panelj.controlpoint)

    end

    return x, r, rj, dmagj, m, nhati
end

"""
    get_relative_geometry_axisym(panel, field_point)

Obtain relevant geometry associated with ring singularity influence calculations for arbitrary field point

**Arguments:**
- `panel::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).
- `field_point::Array{Float}` : [x;r] coordinates of the field point in question.

**Returns:**
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `dmagj::Float` : length of the jth panel
- `m::Float` : Elliptic Function parameter
"""
function get_relative_geometry_axisym(panel, field_point)

    #rename for convenience
    dmagj = panel.length

    #panel control point
    xj = panel.controlpoint[1]
    rj = panel.controlpoint[2]

    #field point
    xi = field_point[1]
    ri = field_point[2]

    #get x and r for these panels
    x = (xi - xj) / rj
    r = ri / rj

    #get phi for these panels
    m = 4.0 * r / (x^2 + (r + 1.0)^2)

    return x, r, rj, dmagj, m
end
