# - If single airfoil, need to put Panel object in a vector - #
function generate_system_geometry(method::Lewis, panel_geometry)
    return generate_system_geometry(method, [panel_geometry])
end

function generate_system_geometry(method::Lewis, panel_geometry::AbstractVector)

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
    # Panel Length (contained in panel_geometry objects)
    panel_length = zeros(TF, (total_panels))

    # x-component of normalized distance from influencing panel center to field point
    x = zeros(TF, (total_panels, total_panels))

    # r-component of normalized distance from influencing panel center to field point
    r = zeros(TF, (total_panels, total_panels))

    # variable used in elliptic function calculations
    k2 = zeros(TF, (total_panels, total_panels))

    system_geometry = (; nbodies, panel_indices, mesh2panel, x, r, k2)

    return generate_system_geometry!(method, system_geometry, panel_geometry)
end

function generate_system_geometry!(method::Lewis, system_geometry, panel_geometry)

    # extract fields
    (; nbodies, panel_indices, mesh2panel, x, r, k2) = system_geometry

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panel_geometry --- ###
            for i in panel_indices[m]
                for j in panel_indices[n]

                    # Get x-locations of influencing and influenced panel_geometry
                    xi = panel_geometry[m].panel_center[mesh2panel[i], 1]
                    xj = panel_geometry[n].panel_center[mesh2panel[j], 1]

                    # Get r-locations of influencing and influenced panel_geometry
                    ri = panel_geometry[m].panel_center[mesh2panel[i], 2]
                    rj = panel_geometry[n].panel_center[mesh2panel[j], 2]

                    # Calculate normalized distance components for current set of panel_geometry
                    x[i, j] = (xi - xj) / rj
                    r[i, j] = ri / rj

                    # Calculate the k^2 value for the elliptic integrals
                    k2[i, j] = 4.0 * r[i, j] / (x[i, j]^2 + (r[i, j] + 1.0)^2)
                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    return system_geometry
end
