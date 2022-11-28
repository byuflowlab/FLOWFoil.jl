#=

Meshing Functions

Authors: Judd Mehr,

=#

abstract type Mesh end

######################################################################
#                                                                    #
#                           PLANAR MESHES                            #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    PlanarMesh{TF,TB,TN}

Mesh for single body.

**Fields:**
 - `nodes::Array{Array{Float,2}}` : [x y] node (panel edge) locations for airfoil
 - `chord::Float` : airfoil chord length
 - `blunt_te::Bool` : boolean for whether or not the trailing edge is blunt or not.
 - `trailing_edge_gap::Float` : trailing edge gap distance
 - `tdp::Float` : dot product of unit vectors of trailing edge bisection and gap vectors
 - `txp::Float` : pseudo-cross product of unit vectors of trailing edge bisection and gap vectors
**Assuptions:**
 - x and y coordinates start at the bottom trailing edge and proceed clockwise.

"""
struct PlanarMesh{TF,TB,TN<:Vector{Matrix{TF}}} <: Mesh
    nodes::TN
    chord::TF
    blunt_te::TB
    trailing_edge_gap::TF
    tdp::TF
    txp::TF
end

"""
    PlanarMeshSystem{TM,TF,TL}

System of meshes to solve.

**Fields:**
 - `meshes::Array{Mesh}` : Array of mesh objects.
 - `scales::Vector{Float}` : Airfoil scaling factors.
 - `angles::Vector{Float}` : Airfoil angles of attack.
 - `locations::Array{Array{TF}}` : Array of leading edge locations.

"""
struct PlanarMeshSystem{TM,TF,TL<:Vector{Matrix{TF}}}
    meshes::TM
    scales::TF
    angles::TF
    locations::TL
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

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

######################################################################
#                                                                    #
#                       AXISYMMETRIC MESHES                          #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    AxiSymMesh{TP,TB}

Axisymmetric Mesh Object

**Fields:**
- `panels::FLOWFoil.AxiSymPanel` : panel objects describing surface geometry.
- `bodyofrevolution::Bool` : Flag as to whether or not the mesh represents a body of revolution.
"""
struct AxiSymMesh{TP,TB} <: Mesh
    panels::TP
    bodyofrevolution::TB
end

"""
    AxiSymPanel{TF,TA}

Panel object for axisymmetric meshes.

**Fields:**
- `controlpoint::Array{Float}` : [x;r] coordinates of panel midpoint.
- `length::Float` : length of panel
- `normal::Array{Float}` : unit normal vector of panel (TODO: remove if unused)
- `beta::Float` : angle panel makes with positive x-axis (radians)
- `radiusofcurvature::Float` : the radius of curvature of the geometry at the panel control point. TODO: make sure this is actually correct with current implementation.
"""
struct AxiSymPanel{TF,TA}
    controlpoint::TA
    length::TF
    normal::TA
    beta::TF
    radiusofcurvature::TF
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

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

######################################################################
#                                                                    #
#                          PERIODIC MESHES                           #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#
