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

