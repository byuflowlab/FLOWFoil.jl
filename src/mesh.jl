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

# TODO: There are probably a bunch of places in this function that broadcasting could be leveraged to shorten code/simplify things.
function generate_mesh(p::PlanarProblem, panels; gap_tolerance=1e-10)

    ### --- Convenience Variables --- ###
    nbodies = length(panels)
    npanels = [panels[i].npanels for i in 1:nbodies]
    total_panels = sum(npanels)
    nnodes = [panels[i].npanels + 1 for i in 1:nbodies]
    total_nodes = sum(nnodes)

    # - Define Body Indexing - #

    #find starting indices for each body
    csnodes = cumsum(nnodes)
    cspanels = cumsum(npanels)

    # put together index ranges of panels for each body
    node_indices = [
        (1 + (i == 1 ? 0 : csnodes[i - 1])):(csnodes[i]) for i in 1:length(nbodies)
    ]
    panel_indices = [
        (1 + (i == 1 ? 0 : cspanels[i - 1])):(cspanels[i]) for i in 1:length(nbodies)
    ]

    ### --- Initialize Vectors --- ###
    TF = typeof(sum([panels[i].panel_length[1] for i in 1:nbodies]))

    # Individual Chord Lengths
    leading_edges = zeros(TF, nbodies)
    trailing_edges = zeros(TF, nbodies)

    ### --- General Mesh Fields --- ###
    # Panel Length (contained in panels objects)
    panel_length = zeros(TF, (total_panels))

    # Distance from influencing panel edge 1 to field points (panel edges)
    r1 = zeros(TF, (total_nodes, total_panels))
    # # Distance from influencing panel edge 2 to field points (panel edges)
    r2 = zeros(TF, (total_nodes, total_panels))

    # Portion of distance 1 normal to the panel of influence
    r1normal = zeros(TF, (total_nodes, total_panels))
    # Portion of distance 1 tangent to the panel of influence
    r1tangent = zeros(TF, (total_nodes, total_panels))

    # Natural log of distance 1
    lnr1 = zeros(TF, (total_nodes, total_panels))
    # # Natural log of distance 2
    lnr2 = zeros(TF, (total_nodes, total_panels))

    # Angle between panel and field point with vertex at edge 1 of influencing panel
    theta1 = zeros(TF, (total_nodes, total_panels))
    # # Angle between panel and field point with vertex at edge 2 of influencing panel
    theta2 = zeros(TF, (total_nodes, total_panels))

    ### --- Trailing Edge Mesh Fields --- ###

    # Panel Length (contained in panels objects)
    TE_panel_length = zeros(TF, nbodies)

    # Distance from influencing panel edge 1 to field points (panel edges)
    r1_TE = zeros(TF, (total_nodes, nbodies))
    # # Distance from influencing panel edge 2 to field points (panel edges)
    r2_TE = zeros(TF, (total_nodes, nbodies))

    # Portion of distance 1 normal to the panel of influence
    r1normal_TE = zeros(TF, (total_nodes, nbodies))
    # Portion of distance 1 tangent to the panel of influence
    r1tangent_TE = zeros(TF, (total_nodes, nbodies))

    # Natural log of distance 1
    lnr1_TE = zeros(TF, (total_nodes, nbodies))
    # # Natural log of distance 2
    lnr2_TE = zeros(TF, (total_nodes, nbodies))

    # Angle between panel and field point with vertex at edge 1 of influencing panel
    theta1_TE = zeros(TF, (total_nodes, nbodies))
    # # Angle between panel and field point with vertex at edge 2 of influencing panel
    theta2_TE = zeros(TF, (total_nodes, nbodies))

    # Blunt Trailing Edge Flags
    blunt_te = zeros(Bool, nbodies)

    # tdp values
    tdp = zeros(TF, nbodies)

    # txp values
    txp = zeros(TF, nbodies)

    # trailing edge gap values
    trailing_edge_gap = zeros(TF, nbodies)

    ##### ----- Loop through each of the bodies to be influenced ----- #####
    for m in 1:nbodies

        ### --- Get Body Information --- ###

        # Get trailing edge information for each body
        tdp[m], txp[m], trailing_edge_gap[m], TE_panel_edges, TE_panel_vector, TE_panel_length[m] = get_trailing_edge_info(
            panels[m].panel_edges
        )

        # Get chord length
        leading_edges[m] = minimum(panels[m].panel_edges[:, :, 1])
        trailing_edges[m] = maximum(panels[m].panel_edges[:, :, 1])

        # Check if trailing edge is sharp or not
        if abs(trailing_edge_gap[m]) >
            gap_tolerance * (trailing_edges[m] - leading_edges[m])

            # set blunt_te to true
            blunt_te[m] = true

        else #(closed trailing edge)

            #set blunt_te to false
            blunt_te[m] = false
        end

        # Copy over panel lengths for convenience
        for pi in panel_indices
            panel_length[pi] = panels[m].panel_length
        end

        ##### ----- Loop through each of the bodies doing the influencing ----- #####
        for n in 1:nbodies

            ### --- Loop through each field node (panel edge) in body m --- ###
            for i in node_indices[m]

                # Panel edges of panels being influenced (body m)
                # NOTE: index i goes 1 beyond length of number of panels, so need to repeat over last panel twice
                panidx = i == node_indices[m][end] ? i - 1 : i
                field_panel_edge = panels[m].panel_edges[panidx, :, :]

                # Get vector and magnitude from first edge of the panel of influence to the field point (edge of panel being influenced)
                # NOTE: index i goes 1 beyond length of number of panels, so need to repeat over last panel twice, using the second panel edge on the repeat
                edgeidx = i == node_indices[m][end] ? 2 : 1

                ### --- Assemble Trailing Edge Gap Panel Infuences --- ###

                # Get trailing edge panel geometry

                # Calculate Influence Geometry
                r1_TE[i, m], r2_TE[i, m], r1normal_TE[i, m], r1tangent_TE[i, m], theta1_TE[i, m], theta2_TE[i, m], lnr1_TE[i, m], lnr2_TE[i, m] = calculate_influence_geometry(
                    TE_panel_edges,
                    TE_panel_vector,
                    TE_panel_length[m],
                    field_panel_edge[edgeidx, :];
                    gap_tolerance=gap_tolerance,
                )

                ### --- Loop through each influence panel in body n --- ###
                for j in panel_indices[n]

                    # - Rename For Convenience - #
                    # Panel edges of influencing panels (body n)
                    influence_panel_edge = panels[n].panel_edges[j, :, :]

                    # Panel Vector and Length
                    influence_panel_vector = panels[n].panel_vector[j, :]

                    # - Calculate Influence Geometry - #

                    r1[i, j], r2[i, j], r1normal[i, j], r1tangent[i, j], theta1[i, j], theta2[i, j], lnr1[i, j], lnr2[i, j] = calculate_influence_geometry(
                        influence_panel_edge,
                        influence_panel_vector,
                        panel_length[j],
                        field_panel_edge[edgeidx, :];
                        gap_tolerance=gap_tolerance,
                    )
                end #for panels being influenced
            end #for influencing panels
        end #for body n
    end #for body m

    # Define System Chord Length
    # TODO: figure out how to expose different options here to the user in a slick way
    chord_length = maximum(trailing_edges) - minimum(leading_edges)

    # - Generate mesh objects - #
    #= TODO:
        There has got to be a better way to do all this.
        The problem is that it's unknown beforehand if there is a trailing edge gap.
        Possibly do the trailing edge gap panel discovery and geometry in the paneling function(s), and then pass those in as part of everything.
        Indexing could get tricky trying to combine everything together though, since the TE panels only influence and are not influenced, so they don't add any more equations to the system.
    =#
    # main mesh
    mesh = PlanarMesh(
        nbodies,
        panel_indices,
        node_indices,
        chord_length,
        panel_length,
        r1,
        lnr1,
        r1normal,
        r1tangent,
        theta1,
        r2,
        lnr2,
        theta2,
    )

    # trailing edge gap panel mesh
    TEmesh = PlanarBluntTEMesh(
        blunt_te,
        trailing_edge_gap,
        tdp,
        txp,
        TE_panel_length,
        r1_TE,
        lnr1_TE,
        r1normal_TE,
        r1tangent_TE,
        theta1_TE,
        r2_TE,
        lnr2_TE,
        theta2_TE,
    )

    return mesh, TEmesh
end

"""
    calculate_influence_geometry(
        influence_panel_edge,
        influence_panel_vector,
        influence_panel_length,
        field_panel_edge;
        gap_tolerance=1e-10,
    )

Calculate all the various geometric pieces of the influence coefficients.

**Arguments:**
- `influence_panel_edge::Array{Float}` : Array of the edge locations of the panel doing the influencing.
- `influence_panel_vector::Vector{Float}` : Vector from the first to second panel edge.
- `influence_panel_length::Float` : Length of the influencing panel.
- `field_point::Vector::Float` : position of field point (point being influence).

**Keyword Arguments:**
- `gap_tolerance::Float` : Distance at which self-induction is assumed (to avoid NaNs and Infs). default = 1e-10

**Returns:**
- `r1::Matrix{TF}` : distance from first panel edge of influence panel to field point
- `r2::Matrix{TF}` : distance from second panel edge of influence panel to field point.
- `r1normal::Matrix{TF}` : component of r1 in normal direction of influence panel
- `r1tangent::Matrix{TF}` : component of r1 in tangent direction of influence panel
- `theta1::Matrix{TF}` : angle between influece panel vector and the r1 vector, adjusted to π for self-influence of panels
- `theta2::Matrix{TF}` angle between influence panel vector and the r2 vector, adjusted to π for self-influence of panels
- `lnr1::Matrix{TF}` : ln(r1), adjusted to zero for self-influence of panels
- `lnr2::Matrix{TF}` : ln(r2), adjusted to zero for self-influence of panels
"""
function calculate_influence_geometry(
    influence_panel_edge,
    influence_panel_vector,
    influence_panel_length,
    field_point;
    gap_tolerance=1e-10,
)

    # Get vector and magnitude from first edge of the panel of influence to the field point (edge of panel being influenced)
    r1vec, r1 = get_r(influence_panel_edge[1, :], field_point)

    # Get vector and magnitude from second edge of the panel of influence to the field point (edge of panel being influenced)
    # NOTE: index i goes 1 beyond length of number of panels, so need to repeat over last panel twice, using the second panel edge on the repeat
    r2vec, r2 = get_r(influence_panel_edge[2, :], field_point)

    # Get the component of r1vec normal to the panel of influence
    r1normal = get_r_normal(r1vec, influence_panel_vector, influence_panel_length)

    # Get the component of r1vec tangent to the panel of influence
    r1tangent = get_r_tangent(r1vec, influence_panel_vector, influence_panel_length)

    # Get the angle between the panel of influence and field point with angle vertex at the first edge of the panel of influence
    theta1 = get_theta(r1normal, r1tangent)

    # Get the angle between the panel of influence and field point with angle vertex at the second edge of the panel of influence
    theta2 = get_theta2(r1normal, r1tangent, influence_panel_length)

    # Calculate natural log values
    #= Natural Log values are to be used in calculating the influence coefficients
        We set the natural logs to zero if the field point is very close to the panel,
        that is, for the self-influence cases.
        NOTE: probably will have other probems if another panel is as close as a panel is to itself...
        We also adjust the values of the angles previously calculated for self-induction cases to be pi or zero depending on positions.
    =#
    if r1 < gap_tolerance
        lnr1 = 0.0
        theta1 = pi
        theta2 = pi

    else
        lnr1 = log(r1)
    end #if

    if r2 < gap_tolerance
        lnr2 = 0.0
        theta1 = 0.0
        theta2 = 0.0

    else
        lnr2 = log(r2)
    end #if

    return r1, r2, r1normal, r1tangent, theta1, theta2, lnr1, lnr2
end

# - If single airfoil, need to put Panel object in a vector - #
function generate_mesh(p::PlanarProblem, panels::TP; gap_tolerance=1e-10) where {TP<:Panel}
    return generate_mesh(p, [panels]; gap_tolerance=gap_tolerance)
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
    r = zeros(TF, (total_panels, total_panels))

    # variable used in elliptic function calculations
    k2 = zeros(TF, (total_panels, total_panels))

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
                    ri = panels[m].panel_center[i, 2]
                    rj = panels[n].panel_center[j, 2]

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
    return AxisymmetricMesh(nbodies, panel_indices, x, r, k2)
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
