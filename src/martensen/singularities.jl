"""
    calculate_periodic_vortex_influence(paneli, panelj)

Cacluate the influence of a periodic vortex at panel j onto panel i.
TODO: BUGGY, NOT WORKING, NEED TO FIX

# Arguments:
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).

# Returns:
- `aij::Float` : Influence of vortex strength at panel j onto panel i.
"""
function calculate_periodic_vortex_influence(
    paneli, panelj, system_geometry, i, j, cascade_parameters
)

    # - Rename for convenience - #
    # control points
    x_i = paneli.panel_center[i]
    y_i = paneli.panel_center[i]

    x_j = panelj.panel_center[j]
    y_j = panelj.panel_center[j]

    # cascade stuff
    stagger = cascade_parameters.stagger
    pitch = cascade_parameters.pitch

    # adjusted angles
    sine_vector_i = sin(stagger + paneli.panel_angle[i])
    cosine_vector_i = cos(stagger + paneli.panel_angle[i])

    sine_vector_j = sin(stagger + paneli.panel_angle[j])
    cosine_vector_j = cos(stagger + paneli.panel_angle[j])

    # - Coefficient - #
    if pitch < eps()
        a = b = k = 0.0
    else
        a = ((x_i - x_j) * cos(stagger) - (y_i - y_j) * sin(stagger)) * 2.0 * pi / pitch
        b = ((x_i - x_j) * sin(stagger) + (y_i - y_j) * cos(stagger)) * 2.0 * pi / pitch
        k = 0.5 / pitch / (cosh(a) - cos(b))
    end

    return (sinh(a) * sine_vector_j - sin(b) * cosine_vector_j) * k * paneli.panel_length[i]
end

function calculate_planar_vortex_influence(paneli, panelj, system_geometry, i, j)

    # rename for convenience
    x_i = paneli.panel_center[i]
    x_j = panelj.panel_center[j]
    y_i = paneli.panel_center[i]
    y_j = panelj.panel_center[j]

    # r = (x_j - x_i)^2 + (y_j - m)^2
    r = system_geometry.x[i, j]^2 + system_geometry.y[i, j]^2

    # u = (y_j - y_i) / (r * 2 * pi)
    u = system_geometry.x[i, j] / (2.0 * pi * r)

    # v = -(x_j - x_i) / (r * 2 * pi)
    v = system_geometry.y[i, j] / (2.0 * pi * r)

    return (u * panelj.cosine_vector[j] + v * panelj.sine_vector[j]) *
           paneli.panel_length[i]
end

function calculate_periodic_self_vortex_influence(panel, i)
    #compute self-inducing coupling coefficients
    # return -0.5 - 2.0 * (panel.delta_angle[i] - 2.0 * pi) / (8 * pi)
    return -0.5 - (panel.delta_angle[i] - pi) / (4 * pi)
end
