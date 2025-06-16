"""
    generate_system_geometry(method::Lewis, panel_geometry)

Generate system geometry for a single airfoil by wrapping it in a vector and
calling the vector version of the function.

# Arguments
- `method::Lewis`: Marker type indicating the Lewis method.
- `panel_geometry`: A single panel geometry object.

# Returns
- `system_geometry`: A NamedTuple containing system-wide geometric information,
  suitable for the Lewis method.
"""
function generate_system_geometry(method::Lewis, panel_geometry)
    return generate_system_geometry(method, [panel_geometry])
end

"""
    generate_system_geometry(method::Lewis, panel_geometry::AbstractVector)

Generate the system geometry for multiple bodies (airfoils) for the Lewis method.

# Arguments
- `method::Lewis`: Marker type indicating the Lewis method.
- `panel_geometry`: Vector of panel geometry objects, each representing one body.

# Returns
- `system_geometry::NamedTuple` with fields:
  - `nbodies::Int`: Number of bodies.
  - `panel_indices::Vector{UnitRange}`: Index ranges for panels of each body.
  - `mesh2panel::Vector{Int}`: Mapping from mesh indices to panel indices.
  - `z::Matrix{TF}`: Normalized axial (z) distances between panel centers.
  - `r::Matrix{TF}`: Normalized radial (r) distances between panel centers.
  - `k2::Matrix{TF}`: k² values used in elliptic integral calculations.
"""
function generate_system_geometry(method::Lewis, panel_geometry::AbstractVector)

    ### --- Convenience Variables --- ###
    nbodies = length(panel_geometry)
    npanels = [pg.npanels for pg in panel_geometry]
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
    # Panel Length (contained in panel_geometry objects)
    panel_length = zeros(TF, (total_panels))

    # z-component of normalized distance from influencing panel center to field point
    z = zeros(TF, (total_panels, total_panels))

    # r-component of normalized distance from influencing panel center to field point
    r = zeros(TF, (total_panels, total_panels))

    # variable used in elliptic function calculations
    k2 = zeros(TF, (total_panels, total_panels))

    system_geometry = (; nbodies, panel_indices, mesh2panel, z, r, k2)

    return generate_system_geometry!(method, system_geometry, panel_geometry)
end

"""
    generate_system_geometry!(method::Lewis, system_geometry, panel_geometry)

Populate and update the `system_geometry` NamedTuple fields by computing the normalized
distances and elliptic integral parameters between panels of all bodies.

# Arguments
- `method::Lewis`: Marker type indicating the Lewis method.
- `system_geometry`: NamedTuple created by `generate_system_geometry`, containing
  pre-allocated arrays and indexing information.
- `panel_geometry`: Vector of panel geometry objects.

# Returns
- The updated `system_geometry` NamedTuple with filled-in distance and elliptic integral parameters.
"""
function generate_system_geometry!(method::Lewis, system_geometry, panel_geometry)

    # extract fields
    (; nbodies, panel_indices, mesh2panel, z, r, k2) = system_geometry

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panel_geometry --- ###
            for i in panel_indices[m]
                for j in panel_indices[n]

                    # Get z-locations of influencing and influenced panel_geometry
                    zi = panel_geometry[m].panel_center[mesh2panel[i], 1]
                    zj = panel_geometry[n].panel_center[mesh2panel[j], 1]

                    # Get r-locations of influencing and influenced panel_geometry
                    ri = panel_geometry[m].panel_center[mesh2panel[i], 2]
                    rj = panel_geometry[n].panel_center[mesh2panel[j], 2]

                    # Calculate normalized distance components for current set of panel_geometry
                    z[i, j] = (zi - zj) / rj
                    r[i, j] = ri / rj

                    # Calculate the k^2 value for the elliptic integrals
                    k2[i, j] = 4.0 * r[i, j] / (z[i, j]^2 + (r[i, j] + 1.0)^2)
                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    return system_geometry
end
