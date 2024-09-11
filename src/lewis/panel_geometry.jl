function generate_panels(p::AxisymmetricProblem, coordinates)

    #broadcast for multiple airfoils
    return generate_panels.(Ref(p), coordinates)
end

function generate_panels(::AxisymmetricProblem, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Separate coordinates
    x = coordinates[:, 1]
    r = coordinates[:, 2]

    # Check if any r coordinates are negative (not allowed in axisymmetric method)
    @assert all(x -> x >= -eps(), r)

    # - Rename for Convenience - #
    npanels = length(x) - 1

    # - Initialize Outputs - #
    panel_center = zeros(TF, npanels, 2)
    panel_length = zeros(TF, npanels)
    panel_normal = zeros(TF, npanels, 2)
    panel_curvature = zeros(TF, npanels)
    panel_angle = zeros(TF, npanels)

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panel_center[i, :] = [0.5 * (x[i] + x[i + 1]); 0.5 * (r[i] + r[i + 1])]

        # Calculate panel length
        panel_vector, panel_length[i] = get_d([x[i] r[i]; x[i + 1] r[i + 1]])

        # Calculate panel unit normal
        panel_normal[i, :] = get_panel_normal(panel_vector, panel_length[i])

        # - Calculate Panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        if panel_vector[1] == 0.0
            beta = pi / 2.0
        else
            beta = atan(panel_vector[2] / panel_vector[1])
        end

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

        #= NOTE:
        #This is the version from the book code.  Maybe it's more robust?

        # sine[i] = (r[i + 1]- r[i]) / panel_length[i]
        sine[i] = (panel_vector[2]) / panel_length[i]
        cosine[i] = (panel_vector[1]) / panel_length[i]
        # cosine[i] = (x[i + 1]- x[i]) / panel_length[i]
        abscos = abs(cosine[i])
        if abscos > ex
            #use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
            beta = atan(sine[i] / cosine[i])
        end

        #if the panel is nearly vertical, set the panel panel_angle to vertical in the correct direction.
        if abscos < ex
            panel_angle[i] = sign(sine[i]) * pi / 2.0
        end

        #otherwise (in most cases)
        if cosine[i] > ex
            panel_angle[i] = beta
        end

        ## For special cases
        #find minimum x point (i.e. the leading edge point
        _, minx = findmin(x)

        #if panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
        if (cosine[i] < -ex) && (i > minx)
            panel_angle[i] = beta - pi
        end

        #if panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
        if (cosine[i] < -ex) && (i < minx)
            panel_angle[i] = beta + pi
        end
        =#

    end

    # - Calculate Panel Curvature - #
    # - If not the end panels, calculate non-zero curvatures - #
    # NOTE: end panels are trailing edges, which are assumed to have zero curvature.
    for i in 2:(npanels - 1)
        panel_curvature[i] = (panel_angle[i + 1] - panel_angle[i - 1]) / 8.0 / pi
    end

    # - Return Panel Object - #
    return AxisymmetricFlatPanel(
        npanels, panel_center, panel_length, panel_normal, panel_angle, panel_curvature
    )
end

function generate_panels!(::AxisymmetricProblem, panels, coordinates::Matrix{TF}) where {TF}

    ### --- SETUP --- ###

    # Separate coordinates
    x = coordinates[:, 1]
    r = coordinates[:, 2]

    # Check if any r coordinates are negative (not allowed in axisymmetric method)
    @assert all(x -> x >= -eps(), r)

    # - Rename for Convenience - #
    npanels = length(x) - 1

    ### --- Loop Through Coordinates --- ###
    for i in 1:npanels

        # Calculate control point (panel center)
        panels.panel_center[i, :] .= [0.5 * (x[i] + x[i + 1]); 0.5 * (r[i] + r[i + 1])]

        # Calculate panels.panel length
        panel_vector, panels.panel_length[i] = get_d([x[i] r[i]; x[i + 1] r[i + 1]])

        # Calculate panels.panel unit normal
        panels.panel_normal[i, :] = get_panel_normal(panel_vector, panels.panel_length[i])

        # - Calculate panels.panel Angles - #
        # Find minimum x point (i.e. the leading edge point) to distinguish between top and bottome of airfoil
        _, minx = findmin(x)

        # NOTE: use standard atan rather than atan2.  For some reason atan2 is not giving the correct angles we want.
        if panel_vector[1] == 0.0
            beta = pi / 2.0
        else
            beta = atan(panel_vector[2] / panel_vector[1])
        end

        # Apply corrections as needed based on orientation of panels.panel in coordinate frame.
        if (panel_vector[1] < 0.0) && (i > minx)
            #if panels.panel is on the top half of the airfoil and has a negative x direction, need to correct the angle from atan
            panels.panel_angle[i] = beta - pi

        elseif (panel_vector[1] < 0.0) && (i < minx)
            #if panels.panel is on the bottom half of the airfoil and has a negative x direction, need to correct the angle from atan
            panels.panel_angle[i] = beta + pi
        else
            panels.panel_angle[i] = beta
        end
    end

    # - Calculate Panel Curvature - #
    # - If not the end panels, calculate non-zero curvatures - #
    # NOTE: end panels are trailing edges, which are assumed to have zero curvature.
    for i in 2:(npanels - 1)
        panels.panel_curvature[i] =
            (panels.panel_angle[i + 1] - panels.panel_angle[i - 1]) / 8.0 / pi
    end

    return nothing
end
