"""
    get_ring_geometry(paneli, panelj)

Obtain relevant geometry associated with ring singularity influence calculations.

# Arguments:
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).

# Returns:
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `dmagj::Float` : length of the jth panel
- `m::Float` : Elliptic Function parameter
- `nhati::Array{Float}` : unit normal vector of ith panel
"""
function get_ring_geometry(paneli, panelj)

    #rename for convenience
    dmagj = panelj.length

    nhati = paneli.normal

    xi = paneli.controlpoint[1]
    ri = paneli.controlpoint[2]

    xj = panelj.controlpoint[1]
    rj = panelj.controlpoint[2]

    #get x and r for these panels
    x = (xi - xj) / rj
    r = ri / rj

    #get phi for these panels
    m = 4.0 * r / (x^2 + (r + 1.0)^2)

    return x, r, rj, dmagj, m, nhati
end

"""
    get_relative_geometry_axisym(panel, field_point)

Obtain relevant geometry associated with ring singularity influence calculations for arbitrary field point

# Arguments:
- `panel::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).
- `field_point::Array{Float}` : [x;r] coordinates of the field point in question.

# Returns:
- `x::Float` : ratio of difference of ith and jth panel x-locations and jth panel r-location ( (xi-xj)/rj )
- `r::Float` : ratio of r-locations of ith and jth panels (ri/rj)
- `rj::Float` : r-location of the jth panel control point
- `dmagj::Float` : length of the jth panel
- `m::Float` : Elliptic Function parameter
"""
function get_relative_geometry_axisym(panel, field_point)

    #rename for convenience
    dmagj = panel.length

    #panel control point
    xj = panel.controlpoint[1]
    rj = panel.controlpoint[2]

    #field point
    xi = field_point[1]
    ri = field_point[2]

    #get x and r for these panels
    x = (xi - xj) / rj
    r = ri / rj

    #get phi for these panels
    m = 4.0 * r / (x^2 + (r + 1.0)^2)

    return x, r, rj, dmagj, m
end
