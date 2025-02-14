function generate_system_geometry(method::Martensen, panels)
    return generate_system_geometry(method, [panels])
end

function generate_system_geometry(method::Martensen, panels::AbstractArray)

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
    return PeriodicMesh(nbodies, panel_indices, x, y, method.pitch, method.stagger)
end
