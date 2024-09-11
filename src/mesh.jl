#=

Meshing Types and Functions

The term Mesh, here, indicates the relational geometry used in calculating influence coefficients for the system.

Authors: Judd Mehr,

=#

######################################################################
#                                                                    #
#                              GENERAL                               #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

abstract type Mesh end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

"""
**Arguments:**
- `p::ProblemType` : Problem type object for dispatch
- `panels::Vector{Panel}` : Array of panel object for airfoil system. (can also be a single panel object if only one body is being modeled)

**Returns:**
- `mesh::Mesh` : Mesh object including various influence geometries for the system
- `TEmesh::Mesh` : Mesh object specifically for trailing edge gap panels if present.
"""
function generate_mesh(p::ProblemType, panels; gap_tolerance=1e-10) end

######################################################################
#                                                                    #
#                           PLANAR MESHES                            #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

#= TODO: This is way too complicated.
The mfoil/xfoil implementation was probably put together for quickness in fortran.
Need to re-derive linear vortex distributions at some point,
  and update to a generalized methodology if possible.
Will probably make the xfoil implementation a special case in its own directory with and Xfoil <: ProblemType dispatch.
=#
"""
    PlanarMesh <: Mesh

Influence geometry for an airfoil or airfoil system.

**Fields:**
- `nbodies::Int64` : number of bodies in the system.
- `panel_indices::Vector{UnitRange{Int64}}` : vector of indices of the overall system matrix associated with each of the panels.
- `node_indices::Vector{UnitRange{Int64}}` : vector of indices of the overall system matrix associated with each panel edge (node).
- `chord::TF` : chord length of system (maximum trailing egde x-coordinate minus minimum leading edge x-coordinate)
- `panel_length::Vector{TF}` : length of each panel
- `r1::Matrix{TF}` : distance from first panel edge of influence panel to field point
- `lnr1::Matrix{TF}` : ln(r1), adjusted to zero for self-influence of panels
- `r1normal::Matrix{TF}` : component of r1 in normal direction of influence panel
- `r1tangent::Matrix{TF}` : component of r1 in tangent direction of influence panel
- `theta1::Matrix{TF}` : angle between influece panel vector and the r1 vector, adjusted to π for self-influence of panels
- `r2::Matrix{TF}` : distance from second panel edge of influence panel to field point.
- `lnr2::Matrix{TF}` : ln(r2), adjusted to zero for self-influence of panels
- `theta2::Matrix{TF}` angle between influence panel vector and the r2 vector, adjusted to π for self-influence of panels

**Assuptions:**
- x and y coordinates start at the bottom trailing edge and proceed clockwise.
"""
struct PlanarMesh{TF} <: Mesh
    nbodies::Int64
    panel_indices::Vector{UnitRange{Int64}}
    node_indices::Vector{UnitRange{Int64}}
    nodes::Matrix{TF}
    mesh2panel::Vector{Int64}
    chord::TF
    panel_length::Vector{TF}
    r1::Matrix{TF}
    lnr1::Matrix{TF}
    r1normal::Matrix{TF}
    r1tangent::Matrix{TF}
    theta1::Matrix{TF}
    r2::Matrix{TF}
    lnr2::Matrix{TF}
    theta2::Matrix{TF}
end

"""
    PlanarBluntTEMesh <: Mesh

Similar to PlanarMesh, but specifically for the trailing edge gap panels (as panels of influence), if they exist.

**Arguments:**
- `blunt_te::Bool` : boolean for whether or not the trailing edge is blunt or not.
- `trailing_edge_gap::Float` : trailing edge gap distance
- `tdp::Float` : dot product of unit vectors of trailing edge bisection and gap vectors
- `txp::Float` : pseudo-cross product of unit vectors of trailing edge bisection and gap vectors
- `panel_length::Vector{TF}` : length of each panel
- `r1::Matrix{TF}` : distance from first panel edge of influence panel to field point
- `lnr1::Matrix{TF}` : ln(r1), adjusted to zero for self-influence of panels
- `r1normal::Matrix{TF}` : component of r1 in normal direction of influence panel
- `r1tangent::Matrix{TF}` : component of r1 in tangent direction of influence panel
- `theta1::Matrix{TF}` : angle between influece panel vector and the r1 vector, adjusted to π for self-influence of panels
- `r2::Matrix{TF}` : distance from second panel edge of influence panel to field point.
- `lnr2::Matrix{TF}` : ln(r2), adjusted to zero for self-influence of panels
- `theta2::Matrix{TF}` angle between influence panel vector and the r2 vector, adjusted to π for self-influence of panels
"""
struct PlanarBluntTEMesh{TF} <: Mesh
    blunt_te::Vector{Bool}
    trailing_edge_gap::Vector{TF}
    tdp::Vector{TF}
    txp::Vector{TF}
    panel_length::Vector{TF}
    r1::Matrix{TF}
    lnr1::Matrix{TF}
    r1normal::Matrix{TF}
    r1tangent::Matrix{TF}
    theta1::Matrix{TF}
    r2::Matrix{TF}
    lnr2::Matrix{TF}
    theta2::Matrix{TF}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

######################################################################
#                                                                    #
#                       AXISYMMETRIC MESHES                          #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    AxisymmetricMesh{TF}

Axisymmetric Mesh Object

**Fields:**
- `nbodies::Int64` : number of bodies in the system.
- `panel_indices::Vector{UnitRange{Int64}}` : vector of indices of the overall system matrix associated with each of the panels.
- `x::Matrix{TF}` : Normalized x-components of distance between panel centers
- `r::Matrix{TF}` : Normalized r-components of distance between panel centers
- `m::Matrix{TF}` : Input values to Elliptic Integral functions
"""
struct AxisymmetricMesh{TF} <: Mesh
    nbodies::Int
    panel_indices::Vector{UnitRange{Int64}}
    mesh2panel::Vector{Int64}
    x::Matrix{TF}
    r::Matrix{TF}
    m::Matrix{TF}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

# - If single airfoil, need to put Panel object in a vector - #
function generate_mesh(p::AxisymmetricProblem, panels::TP; ex=1e-5) where {TP<:Panel}
    return generate_mesh(p, [panels]; ex=ex)
end

function generate_mesh(axisym::AxisymmetricProblem, panels; ex=1e-5)

    ### --- Convenience Variables --- ###
    nbodies = length(panels)
    npanels = [panels[i].npanels for i in 1:nbodies]
    total_panels = sum(npanels)

    # - Define Body Indexing - #

    #find starting indices for each body
    cspanels = cumsum(npanels)

    # put together index ranges of panels for each body
    panel_indices = [
        (1 + (i == 1 ? 0 : cspanels[i - 1])):(cspanels[i]) for i in 1:nbodies
    ]

    # - Map indices - #
    mesh2panel = reduce(vcat,[1:npanels[i] for i in 1:nbodies])

    ### --- Initialize Vectors --- ###
    TF = typeof(sum([panels[i].panel_length[1] for i in 1:nbodies]))

    ### --- General Mesh Fields --- ###
    # Panel Length (contained in panels objects)
    panel_length = zeros(TF, (total_panels))

    # x-component of normalized distance from influencing panel center to field point
    x = zeros(TF, (total_panels, total_panels))

    # r-component of normalized distance from influencing panel center to field point
    r = zeros(TF, (total_panels, total_panels))

    # variable used in elliptic function calculations
    k2 = zeros(TF, (total_panels, total_panels))

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panels --- ###
            for i in panel_indices[m]
                for j in panel_indices[n]

                    # Get x-locations of influencing and influenced panels
                    xi = panels[m].panel_center[mesh2panel[i], 1]
                    xj = panels[n].panel_center[mesh2panel[j], 1]

                    # Get r-locations of influencing and influenced panels
                    ri = panels[m].panel_center[mesh2panel[i], 2]
                    rj = panels[n].panel_center[mesh2panel[j], 2]

                    # Calculate normalized distance components for current set of panels
                    x[i, j] = (xi - xj) / rj
                    r[i, j] = ri / rj

                    # Calculate the k^2 value for the elliptic integrals
                    k2[i, j] = 4.0 * r[i, j] / (x[i, j]^2 + (r[i, j] + 1.0)^2)
                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    return AxisymmetricMesh(nbodies, panel_indices, mesh2panel,x, r, k2)
end

######################################################################
#                                                                    #
#                          PERIODIC MESHES                           #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    PeriodicMesh{TF}

Periodic Mesh Object

**Fields:**
- `nbodies::Int64` : number of bodies in the system.
- `panel_indices::Vector{UnitRange{Int64}}` : vector of indices of the overall system matrix associated with each of the panels.
- `x::Matrix{TF}` : x-components of distance between panel centers
- `y::Matrix{TF}` : y-components of distance between panel centers
- `pitch::Vector{Float}` : Distance between airfoils in cascade
"""
struct PeriodicMesh{TF} <: Mesh
    nbodies::Int
    panel_indices::Vector{UnitRange{Int64}}
    x::Matrix{TF}
    y::Matrix{TF}
    pitch::TF
    stagger::TF
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

# - If single airfoil, need to put Panel object in a vector - #
function generate_mesh(pp::PeriodicProblem, panels::TP) where {TP<:Panel}
    return generate_mesh(pp, [panels])
end

function generate_mesh(pp::PeriodicProblem, panels)

    ### --- Convenience Variables --- ###
    nbodies = length(panels)
    npanels = [panels[i].npanels for i in 1:nbodies]
    total_panels = sum(npanels)

    # - Define Body Indexing - #

    #find starting indices for each body
    cspanels = cumsum(npanels)

    # put together index ranges of panels for each body
    panel_indices = [
        (1 + (i == 1 ? 0 : cspanels[i - 1])):(cspanels[i]) for i in 1:length(nbodies)
    ]

    ### --- Initialize Vectors --- ###
    TF = typeof(sum([panels[i].panel_length[1] for i in 1:nbodies]))

    ### --- General Mesh Fields --- ###
    # Panel Length (contained in panels objects)
    panel_length = zeros(TF, (total_panels))

    # x-component of normalized distance from influencing panel center to field point
    x = zeros(TF, (total_panels, total_panels))

    # r-component of normalized distance from influencing panel center to field point
    y = zeros(TF, (total_panels, total_panels))

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panels --- ###
            for i in panel_indices[m]
                for j in panel_indices[m]

                    # Get x-locations of influencing and influenced panels
                    xi = panels[m].panel_center[i, 1]
                    xj = panels[n].panel_center[j, 1]

                    # Get r-locations of influencing and influenced panels
                    yi = panels[m].panel_center[i, 2]
                    yj = panels[n].panel_center[j, 2]

                    # Calculate normalized distance components for current set of panels
                    x[i, j] = xi - xj
                    y[i, j] = yi - yj
                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    return PeriodicMesh(nbodies, panel_indices, x, y, pp.pitch, pp.stagger)
end
