function generate_system_geometry(method::HessSmith, panel_geometry)
    return generate_system_geometry(method, [panel_geometry])
end

function generate_system_geometry(method::HessSmith, panel_geometry::AbstractVector)

    ### --- Convenience Variables --- ###
    nbodies = length(panel_geometry)
    npanels = [panel_geometry[i].npanels for i in 1:nbodies]
    total_panels = sum(npanels)

    # - Define Body Indexing - #

    # find starting indices for each body
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

    # needed distance
    r_influence = zeros(TF, (total_panels, total_panels))

    system_geometry = (;
        nbodies,
        panel_indices,
        mesh2panel,
        r_x,
        r_y,
        r_squared,
        r_influence
    )

    return generate_system_geometry!(method, system_geometry, panel_geometry)
end

function generate_system_geometry!(
    method::HessSmith, system_geometry, panel_geometry::AbstractVector
)
    (; nbodies, panel_indices, mesh2panel, r_x, r_y, r_squared, r_influence) = system_geometry

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panel_geometry --- ###
            for i in panel_indices[m]
                for j in panel_indices[m]

                    # Get x-locations of influencing and influenced panel_geometry
                    xi_control = panel_geometry[m].panel_center[mesh2panel[i], 1]
                    xj_field = panel_geometry[m].x[mesh2panel[j], 1]

                    # Get y-locations of influencing and influenced panel_geometry
                    yi_control = panel_geometry[m].panel_center[mesh2panel[i], 2]
                    yj_field = panel_geometry[m].y[mesh2panel[j], 2]

                    # Calculate distance components for current set of panel_geometry
                    r_x[i, j] = xj_field - xi_control
                    r_y[i, j] = yj_field - yi_control
                    r_squared[i, j] = r_x[i, j]^2 + r_y[i, j]^2
                    r_influence[i, j] = sqrt(r_squared[i, j])
                end #for jth influencing panel
            end #for ith influenced panel
        end #for nth influencing body
    end #for mth influenced body

    # Return Mesh
    return system_geometry
end
