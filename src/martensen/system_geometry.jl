"""
    generate_system_geometry(method::Martensen, panel_geometry)

Generates system geometry for a single body by wrapping the vector version.

# Arguments
- `method::Martensen`: The panel method object containing method parameters.
- `panel_geometry`: Panel geometry object for a single body.

# Returns
- A system geometry named tuple containing indexing, distance arrays, and pitch.
"""
function generate_system_geometry(method::Martensen, panel_geometry)
    return generate_system_geometry(method, [panel_geometry])
end

"""
    generate_system_geometry(method::Martensen, panel_geometry::AbstractVector)

Generates system geometry data for one or more bodies in a panel method simulation.

# Arguments
- `method::Martensen`: The panel method object containing parameters like solidity.
- `panel_geometry::AbstractVector`: Vector of panel geometry objects, one per body.

# Returns
- A named tuple containing:
  - `nbodies`: Number of bodies.
  - `panel_indices`: Vector of ranges indexing each body's panels within the global panel mesh.
  - `mesh2panel`: Vector mapping mesh indices to panel indices.
  - `r_x`, `r_y`: Matrices containing x and y distances between panel centers.
  - `r_squared`: Matrix of squared distances.
  - `pitch`: Characteristic pitch length for cascade geometry.
"""
function generate_system_geometry(method::Martensen, panel_geometry::AbstractVector)

    ### --- Convenience Variables --- ###
    nbodies = length(panel_geometry)
    npanels = [panel_geometry[i].npanels for i in 1:nbodies]
    total_panels = sum(npanels)

    # - Define Body Indexing - #

    #find starting indices for each body
    cspanels = cumsum(npanels)

    # put together index ranges of panel_geometry for each body
    panel_indices = [(1 + (i == 1 ? 0 : cspanels[i - 1])):(cspanels[i]) for i in 1:nbodies]

    # - Map indices - #
    mesh2panel = reduce(vcat, [1:npanels[i] for i in 1:nbodies])

    ### --- Initialize Vectors --- ###
    TF = typeof(sum([panel_geometry[i].panel_length[1] for i in 1:nbodies]))

    ### --- General Mesh Fields --- ###
    # x-component of distance from influencing panel center to field point
    r_x = zeros(TF, (total_panels, total_panels))

    # y-component of distance from influencing panel center to field point
    r_y = zeros(TF, (total_panels, total_panels))

    # total distance
    r_squared = zeros(TF, (total_panels, total_panels))

    system_geometry = (;
        nbodies,
        panel_indices,
        mesh2panel,
        r_x,
        r_y,
        r_squared,
        pitch=calculate_chord(panel_geometry) / method.solidity,
    )
    return generate_system_geometry!(method, system_geometry, panel_geometry)
end

"""
    generate_system_geometry!(method::Martensen, system_geometry, panel_geometry::AbstractVector)

Populates distance matrices (`r_x`, `r_y`, `r_squared`) within the provided `system_geometry`
based on panel center coordinates for all bodies.

# Arguments
- `method::Martensen`: The panel method object.
- `system_geometry`: Named tuple with pre-allocated fields for geometry data.
- `panel_geometry::AbstractVector`: Vector of panel geometry objects for each body.

# Returns
- The updated `system_geometry` with computed distances filled in.
"""
function generate_system_geometry!(
    method::Martensen, system_geometry, panel_geometry::AbstractVector
)
    (; nbodies, panel_indices, mesh2panel, r_x, r_y, r_squared) = system_geometry

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panel_geometry --- ###
            for i in panel_indices[m]
                for j in panel_indices[m]

                    # Get x-locations of influencing and influenced panel_geometry
                    xi = panel_geometry[m].panel_center[mesh2panel[i], 1]
                    xj = panel_geometry[m].panel_center[mesh2panel[j], 1]

                    # Get y-locations of influencing and influenced panel_geometry
                    yi = panel_geometry[m].panel_center[mesh2panel[i], 2]
                    yj = panel_geometry[m].panel_center[mesh2panel[j], 2]

                    # Calculate distance components for current set of panel_geometry
                    r_x[i, j] = xj - xi
                    r_y[i, j] = yj - yi
                    r_squared[i, j] = r_x[i, j]^2 + r_y[i, j]^2
                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    return system_geometry
end
