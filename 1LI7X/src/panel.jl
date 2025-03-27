#=

Panel Generation Types and Functions

Panels are defined as the geometry of the discretized pieces of the overall geometry, organized in a way that will be helpful for finding geometric relations later (meshing).

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                              GENERAL                               #
#                                                                    #
######################################################################

#---------------------------------#
#              Types              #
#---------------------------------#

abstract type Panel end

#---------------------------------#
#            Functions            #
#---------------------------------#

"""
    generate_panels(p, coordinates)

Generate panel object for a give set of coordinates.

**Arguments:**
- `p::ProblemType` : problem type object
- `coordinates::NTuple{Matrix{Float}}` : Tuple of [x y] matrices of airfoil coordinates (may be a single matrix as well)

**Returns:**
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
"""
function generate_panels(::ProblemType, coordinates) end

######################################################################
#                                                                    #
#                              PLANAR                                #
#                                                                    #
######################################################################

#---------------------------------#
#              Types              #
#---------------------------------#

"""
    LinearFlatPanel

**Fields:**
- `npanels::Vector{Int}` : number of panels on each body
- `panel_edges::Array{Float, 3}` : Panel edge [x z] coordinates, indexed as [panel number, edge (1 or 2), coordinate (x or z)]
- `panel_vector::Matrix{TF}` : Vectors from panel edge 1 to panel edge 2, index as [panel number, edge number]
- `panel_length::Vector{TF}` : Lengths of panels (magnitude of panel vector)
"""
struct LinearFlatPanel{TF} <: Panel
    npanels::Int
    panel_edges::Array{TF,3}
    panel_vector::Matrix{TF}
    panel_length::Vector{TF}
    node::Matrix{TF}
end

#---------------------------------#
#            Functions            #
#---------------------------------#

#= NOTE:
The current implentation here is the Xfoil-like implementation.  Eventually, this will likely be moved to it's own dispatch method as better methods are developed.
=#
function generate_panels(p::PlanarProblem, coordinates)

    #broadcast for multiple airfoils
    return reduce(vcat, generate_panels.(Ref(p), coordinates))
end

function generate_panels(::PlanarProblem, coordinates::Matrix{TF}) where {TF}

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Get number of airfoil nodes for convenience
    numpanels = length(x) - 1

    # Initialize Outputs
    panel_edges = zeros(TF, numpanels, 2, 2)
    panel_lengths = zeros(TF, numpanels)
    panel_vectors = zeros(TF, numpanels, 2)
    nodes = zeros(TF, numpanels + 1, 2)
    nodes[1, :] = [x[1] y[1]]
    nodes[end, :] = [x[end] y[end]]

    for i in 1:numpanels
        # Get node locations from x,y coordinates
        panel_edges[i, :, :] = [x[i] y[i]; x[i + 1] y[i + 1]]

        # Calculate Panel Lengths
        panel_lengths[i] = sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)

        # Calculate Panel Vectors
        panel_vectors[i, :] = [x[i + 1] - x[i] y[i + 1] - y[i]]

        # Set Node Values
        nodes[i, :] = [x[i] y[i]]
    end

    return [LinearFlatPanel(numpanels, panel_edges, panel_vectors, panel_lengths, nodes)]
end

"""
"""
function generate_panels!(::PlanarProblem, panels, coordinates::Matrix{TF}) where {TF}

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Get number of airfoil nodes for convenience
    numpanels = length(x) - 1

    # Initialize Outputs
    panels.nodes[1, :] .= [x[1] y[1]]
    panels.nodes[end, :] .= [x[end] y[end]]

    for i in 1:numpanels
        # Get node locations from x,y coordinates
        panels.panel_edges[i, :, :] .= [x[i] y[i]; x[i + 1] y[i + 1]]

        # Calculate Panel Lengths
        panels.panel_lengths[i] .= sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)

        # Calculate Panel Vectors
        panels.panel_vectors[i, :] .= [x[i + 1] - x[i] y[i + 1] - y[i]]

        # Set Node Values
        panels.nodes[i, :] .= [x[i] y[i]]
    end

    return nothing
end

######################################################################
#                                                                    #
#                            AXISYMMETRIC                            #
#                                                                    #
######################################################################

#---------------------------------#
#              Types              #
#---------------------------------#

"""
    AxisymmetricFlatPanel

**Fields:**
- `npanels::Vector{Int}` : number of panels on each body
- `panel_center::Matrix{TF}` : Panel center locations
- `panel_length::Vector{TF}` : Lengths of panels (magnitude of panel vector)
- `panel_normal::Matrix{TF}` : Panel normal unit vectors
- `panel_angle::Vector{TF}` : Angles of panels
- `panel_curvature::Vector{TF}` : Curvature of panels
"""
struct AxisymmetricFlatPanel{TF} <: Panel
    npanels::Int
    panel_center::Matrix{TF}
    panel_length::Vector{TF}
    panel_normal::Matrix{TF}
    panel_angle::Vector{TF}
    panel_curvature::Vector{TF}
end

#---------------------------------#
#            Functions            #
#---------------------------------#

function generate_panels(p::AxisymmetricProblem, coordinates)

    #broadcast for multiple airfoils
    return generate_panels.(Ref(p), coordinates)
end

function generate_panels(::AxisymmetricProblem, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Separate coordinates
    x = coordinates[:, 1]
    r = coordinates[:, 2]

    # Check if any r coordinates are negative (not allowed in axisymmetric method)
    @assert all(x -> x >= -eps(), r)

    # - Rename for Convenience - #
    npanels = length(x) - 1

    # - Initialize Outputs - #
    panel_center = zeros(TF, npanels, 2)
    panel_length = zeros(TF, npanels)
    panel_normal = zeros(TF, npanels, 2)
    panel_curvature = zeros(TF, npanels)
    panel_angle = zeros(TF, npanels)

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (r[i] + r[i + 1])]

        # Calculate panel length
        panel_vector, panel_length[i] = get_d([x[i] r[i]; x[i + 1] r[i + 1]])

        # Calculate panel unit normal
        panel_normal[i, :] = get_panel_normal(panel_vector, panel_length[i])

        # - Calculate Panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        if panel_vector[1] == 0.0
            beta = pi / 2.0
        else
            beta = atan(panel_vector[2] / panel_vector[1])
        end

        # Apply corrections as needed based on orientation of panel in coordinate frame.
        if (panel_vector[1] < 0.0) && (i > minx)
            #if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta - pi

        elseif (panel_vector[1] < 0.0) && (i < minx)
            #if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta + pi
        else
            panel_angle[i] = beta
        end

        #= NOTE:
        #This is the version from the book code.  Maybe it's more robust?

        # sine[i] = (r[i + 1]- r[i]) / panel_length[i]
        sine[i] = (panel_vector[2]) / panel_length[i]
        cosine[i] = (panel_vector[1]) / panel_length[i]
        # cosine[i] = (x[i + 1]- x[i]) / panel_length[i]
        abscos = abs(cosine[i])
        if abscos > ex
            #use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
            beta = atan(sine[i] / cosine[i])
        end

        #if the panel is nearly vertical, set the panel panel_angle to vertical in the correct direction.
        if abscos < ex
            panel_angle[i] = sign(sine[i]) * pi / 2.0
        end

        #otherwise (in most cases)
        if cosine[i] > ex
            panel_angle[i] = beta
        end

        ## For special cases
        #find minimum x point (i.e. the leading edge point
        _, minx = findmin(x)

        #if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
        if (cosine[i] < -ex) && (i > minx)
            panel_angle[i] = beta - pi
        end

        #if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
        if (cosine[i] < -ex) && (i < minx)
            panel_angle[i] = beta + pi
        end
        =#

    end

    # - Calculate Panel Curvature - #
    # - If not the end panels, calculate non-zero curvatures - #
    # NOTE: end panels are trailing edges, which are assumed to have zero curvature.
    for i in 2:(npanels - 1)
        panel_curvature[i] = (panel_angle[i + 1] - panel_angle[i - 1]) / 8.0 / pi
    end

    # - Return Panel Object - #
    return AxisymmetricFlatPanel(
        npanels, panel_center, panel_length, panel_normal, panel_angle, panel_curvature
    )
end

function generate_panels!(::AxisymmetricProblem, panels, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Separate coordinates
    x = coordinates[:, 1]
    r = coordinates[:, 2]

    # Check if any r coordinates are negative (not allowed in axisymmetric method)
    @assert all(x -> x >= -eps(), r)

    # - Rename for Convenience - #
    npanels = length(x) - 1

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panels.panel_center[i, :] .= [0.5 * (x[i] + x[i + 1]); 0.5 * (r[i] + r[i + 1])]

        # Calculate panels.panel length
        panel_vector, panels.panel_length[i] = get_d([x[i] r[i]; x[i + 1] r[i + 1]])

        # Calculate panels.panel unit normal
        panels.panel_normal[i, :] = get_panel_normal(panel_vector, panels.panel_length[i])

        # - Calculate panels.panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        if panel_vector[1] == 0.0
            beta = pi / 2.0
        else
            beta = atan(panel_vector[2] / panel_vector[1])
        end

        # Apply corrections as needed based on orientation of panels.panel in coordinate frame.
        if (panel_vector[1] < 0.0) && (i > minx)
            #if panels.panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            panels.panel_angle[i] = beta - pi

        elseif (panel_vector[1] < 0.0) && (i < minx)
            #if panels.panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            panels.panel_angle[i] = beta + pi
        else
            panels.panel_angle[i] = beta
        end

    end

    # - Calculate Panel Curvature - #
    # - If not the end panels, calculate non-zero curvatures - #
    # NOTE: end panels are trailing edges, which are assumed to have zero curvature.
    for i in 2:(npanels - 1)
        panels.panel_curvature[i] = (panels.panel_angle[i + 1] - panels.panel_angle[i - 1]) / 8.0 / pi
    end

    return nothing
end

######################################################################
#                                                                    #
#                              PERIODIC                              #
#                                                                    #
######################################################################

#---------------------------------#
#              Types              #
#---------------------------------#

"""
    ConstantFlatPanel

**Fields:**
- `npanels::Vector{Int}` : number of panels on each body
- `panel_center::Matrix{TF}` : Panel center locations
- `panel_length::Vector{TF}` : Lengths of panels (magnitude of panel vector)
- `panel_normal::Matrix{TF}` : Panel normal unit vectors
- `panel_angle::Vector{TF}` : Angles of panels
- `delta_angle::Vector{TF}` : Change in angle of panels from one side to the other
"""
struct ConstantFlatPanel{TF} <: Panel
    npanels::Int
    panel_center::Matrix{TF}
    panel_length::Vector{TF}
    panel_normal::Matrix{TF}
    panel_angle::Vector{TF}
    delta_angle::Vector{TF}
end

#---------------------------------#
#            Functions            #
#---------------------------------#

function generate_panels(p::PeriodicProblem, coordinates)

    #broadcast for multiple airfoils
    return generate_panels.(Ref(p), coordinates)
end

function generate_panels(::PeriodicProblem, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # - Rename for Convenience - #
    npanels = length(x) - 1

    # - Initialize Outputs - #
    panel_center = zeros(TF, npanels, 2)
    panel_length = zeros(TF, npanels)
    panel_normal = zeros(TF, npanels, 2)
    panel_curvature = zeros(TF, npanels)
    panel_angle = zeros(TF, npanels)
    delta_angle = zeros(TF, npanels)

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (y[i] + y[i + 1])]

        # Calculate panel length
        panel_vector, panel_length[i] = get_d([x[i] y[i]; x[i + 1] y[i + 1]])

        # Calculate panel unit normal
        panel_normal[i, :] = get_panel_normal(panel_vector, panel_length[i])

        # - Calculate Panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        beta = atan(panel_vector[2] / panel_vector[1])

        # Apply corrections as needed based on orientation of panel in coordinate frame.
        if (panel_vector[1] < 0.0) && (i > minx)
            #if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta - pi

        elseif (panel_vector[1] < 0.0) && (i < minx)
            #if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta + pi
        else
            panel_angle[i] = beta
        end
    end

    for i in 1:npanels
        if i == 1 || i == npanels
            delta_angle[i] = (panel_angle[i]) / 2.0
        else
            delta_angle[i] = (panel_angle[i + 1] - panel_angle[i - 1]) / 2.0
        end
    end

    # - Return Panel Object - #
    return ConstantFlatPanel(
        npanels, panel_center, panel_length, panel_normal, panel_angle, delta_angle
    )
end
