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

#= TODO: This is way too complicated.
The mfoil/xfoil implementation was probably put together for quickness in fortran.
Need to re-derive linear vortex distributions at some point
  and update to a generalized methodology if possible.
=#
"""
    PlanarMesh{TF,TB,TN}

Mesh for single body.

**Fields:**
- `chord::Float` : airfoil chord length
- `blunt_te::Bool` : boolean for whether or not the trailing edge is blunt or not.
- `trailing_edge_gap::Float` : trailing edge gap distance
- `tdp::Float` : dot product of unit vectors of trailing edge bisection and gap vectors
- `txp::Float` : pseudo-cross product of unit vectors of trailing edge bisection and gap vectors

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

"""
**Arguments:**
- `coordinates::NTuple{N,Matrix{Float}}` : Tuple containting arrays of both x and y coordinates (x first column, y second column) for each airfoil in the airfoil system.
"""
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
    # (see absolute_chord and sum_chord functions below)
    chord_length = maximum(trailing_edges) - minimum(leading_edges)

    # - Generate mesh objects - #
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

function calculate_influence_geometry(
    influence_panel_edge,
    influence_panel_vector,
    influence_panel_length,
    field_panel_edge;
    gap_tolerance=1e-10,
)
    r1vec, r1 = get_r(influence_panel_edge[1, :], field_panel_edge)

    # Get vector and magnitude from second edge of the panel of influence to the field point (edge of panel being influenced)
    # NOTE: index i goes 1 beyond length of number of panels, so need to repeat over last panel twice, using the second panel edge on the repeat
    r2vec, r2 = get_r(influence_panel_edge[2, :], field_panel_edge)

    # Get the component of r1vec normal to the panel of influence
    r1normal = get_r_normal(r1vec, influence_panel_vector, influence_panel_length)

    # Get the component of r1vec tangent to the panel of influence
    r1tangent = get_r_tangent(r1vec, influence_panel_vector, influence_panel_length)

    # Get the angle between the panel of influence and field point with angle vertex at the first edge of the panel of influence
    theta1 = get_theta(r1normal, r1tangent)

    # Get the angle between the panel of influence and field point with angle vertex at the second edge of the panel of influence
    theta2 = get_theta2(r1normal, r1tangent, influence_panel_length)

    # Calculate natural log values to be used in calculating the influence coefficients, setting to zero if self-induction is likely.
    # NOTE: probably will have other probems if another panel is as close as a panel is to itself...
    # Also adjust the values of the angles for self-induction cases.
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

function generate_mesh(p::PlanarProblem, panels::TP; gap_tolerance=1e-10) where {TP<:Panel}
    return generate_mesh(p, [panels]; gap_tolerance=gap_tolerance)
end

"""
"""
function absolute_chord(leading_edges, trailing_edges)
    return maximum(trailing_edges) - minimum(leading_edges)
end

"""
"""
function sum_chord(leading_edges, trailing_edges)
    return sum(trailing_edges .- leading_edges)
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
- `controlpoint::Array{Float}` : [x;r] coordinates of panel midpoint.
- `length::Float` : length of panel
- `normal::Array{Float}` : unit normal vector of panel (TODO: remove if unused)
- `beta::Float` : angle panel makes with positive x-axis (radians)
- `radius_of_curvature::Float` : the radius of curvature of the geometry at the panel control point. TODO: make sure this is actually correct with current implementation.
- `body_of_revolution::Bool` : Flag as to whether or not the mesh represents a body of revolution.
"""
struct AxiSymMesh{TF} <: Mesh
    controlpoint::Vector{Vector{TF}}
    length::Vector{TF}
    normal::Vector{Vector{TF}}
    beta::Vector{TF}
    radius_of_curvature::Vector{TF}
    body_of_revolution::Bool
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

"""
    generate_mesh(x, r; body_of_revolution)

Generate mesh for axisymmetric body.

**Arguments:**
- `coordinates::NTuple{N,Matrix{Float}}` : Tuple containting arrays of both x and y coordinates (x first column, y second column) for each airfoil in the airfoil system.

**Keyword Arguments:**

**Returns:**
- `mesh::FLOWFoil.Array{Mesh}` :
"""
# function generate_mesh(axisym::AxisymmetricProblem, coordinates; ex=1e-5)
# end

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
