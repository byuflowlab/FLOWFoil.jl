function generate_panels(p::PeriodicProblem, coordinates)

    #broadcast for multiple airfoils
    return generate_panels.(Ref(p), coordinates)
end

function generate_panels(::PeriodicProblem, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Separate coordinates
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    # - Rename for Convenience - #
    npanels = length(x) - 1

    # - Initialize Outputs - #
    panel_center = zeros(TF, npanels, 2)
    panel_length = zeros(TF, npanels)
    panel_normal = zeros(TF, npanels, 2)
    panel_curvature = zeros(TF, npanels)
    panel_angle = zeros(TF, npanels)
    delta_angle = zeros(TF, npanels)

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (y[i] + y[i + 1])]

        # Calculate panel length
        panel_vector, panel_length[i] = get_d([x[i] y[i]; x[i + 1] y[i + 1]])

        # Calculate panel unit normal
        panel_normal[i, :] = get_panel_normal(panel_vector, panel_length[i])

        # - Calculate Panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        beta = atan(panel_vector[2] / panel_vector[1])

        # Apply corrections as needed based on orientation of panel in coordinate frame.
        if (panel_vector[1] < 0.0) && (i > minx)
            #if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta - pi

        elseif (panel_vector[1] < 0.0) && (i < minx)
            #if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            panel_angle[i] = beta + pi
        else
            panel_angle[i] = beta
        end
    end

    for i in 1:npanels
        if i == 1 || i == npanels
            delta_angle[i] = (panel_angle[i]) / 2.0
        else
            delta_angle[i] = (panel_angle[i + 1] - panel_angle[i - 1]) / 2.0
        end
    end

    # - Return Panel Object - #
    return ConstantFlatPanel(
        npanels, panel_center, panel_length, panel_normal, panel_angle, delta_angle
    )
end
