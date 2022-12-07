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
    PlanarFlatPanel

**Fields:**
- `npanels::Vector{Int}` : number of panels on each body
- `panel_edges::Array{Float, 3}` : Panel edge [x z] coordinates, indexed as [panel number, edge (1 or 2), coordinate (x or z)]
- `panel_vector::Matrix{TF}` : Vectors from panel edge 1 to panel edge 2, index as [panel number, edge number]
- `panel_length::Vector{TF}` : Lengths of panels (magnitude of panel vector)
"""
struct PlanarFlatPanel{TF} <: Panel
    npanels::Int
    panel_edges::Array{TF,3}
    panel_vector::Matrix{TF}
    panel_length::Vector{TF}
end

#---------------------------------#
#            Functions            #
#---------------------------------#

#= NOTE:
The current implentation here is the Xfoil-like implementation.  Eventually, this will likely be moved to it's own dispatch method as better methods are developed.
=#
function generate_panels(p::PlanarProblem, coordinates)

    #broadcast for multiple airfoils
    return generate_panels.(Ref(p), coordinates)
end

function generate_panels(::PlanarProblem, coordinates::Matrix{TF}) where {TF}

    # Separate out coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # Get number of airfoil nodes for convenience
    numpanels = length(x) - 1

    # Initialize Outputs
    panel_edges = zeros(TF, numpanels, 2, 2)
    panel_lengths = zeros(TF, numpanels)
    panel_vectors = zeros(TF, numpanels, 2)

    for i in 1:numpanels
        # Get node locations from x,y coordinates
        panel_edges[i, :, :] = [x[i] y[i]; x[i + 1] y[i + 1]]

        # Calculate Panel Lengths
        panel_lengths[i] = sqrt((x[i + 1] - x[i])^2 + (y[i + 1] - y[i])^2)

        # Calculate Panel Vectors
        panel_vectors[i, :] = [x[i + 1] - x[i] y[i + 1] - y[i]]
    end
    return PlanarFlatPanel(numpanels, panel_edges, panel_vectors, panel_lengths)
end

######################################################################
#                                                                    #
#                            AXISYMMETRIC                            #
#                                                                    #
######################################################################

#TODO: look at planar struct for updated formatting
"""
"""
struct AxisymmetricFlatPanel{TF} <: Panel
    npanels::Int
    control_point::Vector{Vector{TF}}
    panel_length::Vector{TF}
    panel_normal::Vector{Vector{TF}}
    panel_curvature::Vector{TF}
    panel_angle::Vector{TF}
end

function generate_panels(p::AxisymmetricProblem, coordinates)

    #broadcast for multiple airfoils
    return generate_panels.(Ref(p), coordinates)
end

function generate_panels(::AxisymmetricProblem, coordinates::Matrix{TF}) where {TF}
    # Separate out coordinates
    x = coordinates[:, 1]
    r = coordinates[:, 2]

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
        ## sine[i] = (r[i + 1]- r[i]) / dmag[i]
        #sine[i] = (d[2]) / dmag[i]
        #cosine[i] = (d[1]) / dmag[i]
        ## cosine[i] = (x[i + 1]- x[i]) / dmag[i]
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
        #    slope[i] = beta- pi
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
end
