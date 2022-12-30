#=
Geometry Engine for FLOWFoil.jl

Authors: Judd Mehr,

=#

"""
    get_trailing_edge_info(nodes)

Calculate various items needed for trailing edge treatment.

**Arguments:**
 - `nodes::Array{Float,2}` : Array of [x y] locations for the airfoil nodes.

**Returns:**
 - `tdp::Float` : dot product of TE bisection and TE gap unit vectors
 - `txp::Float` : "cross product" of TE bisection and TE gap unit vectors
 - `trailing_edge_gap::Float` : TE gap distance
"""
function get_trailing_edge_info(panel_edges)

    # - get bisection vector - #
    # get vector along first panel
    d1, d1mag = get_d(panel_edges[1, 2, :], panel_edges[1, 1, :])
    d1 /= d1mag

    # get vector along second panel
    dn, dnmag = get_d(panel_edges[end, :, :])
    dn /= dnmag

    # calculate vector that bisects the first and last panel vectors
    bisector = 0.5 * (d1 + dn)

    # - Calculate tdp - #
    # normalize to get the unit vector
    ttehat = bisector / sqrt(bisector[1]^2 + bisector[2]^2)

    # gap edges
    gap_edges = [panel_edges[end, 2, :]'; panel_edges[1, 1, :]']

    # get panel vector
    dte, dtemag = get_d(panel_edges[1, 1, :], panel_edges[end, 2, :])

    if dtemag == 0.0
        return 1.0, 0.0, 0.0, gap_edges, dte, dtemag
    else
        # normalize panelvector
        dtehat = dte / dtemag

        # get dot product of bisection vector and panel vector.
        tdp = ttehat[1] * dtehat[1] + ttehat[2] * dtehat[2]

        # - Calculate txp - #
        # get cross product of bisection vector and panel vector
        txp = abs(ttehat[1] * dtehat[2] - ttehat[2] * dtehat[1])

        # - Get trailing edge gap - #
        trailing_edge_gap = -dte[1] * ttehat[2] + dte[2] * ttehat[1]

        return tdp, txp, trailing_edge_gap, gap_edges, dte, dtemag
    end
end

"""
    function get_r(node,point)

Calculate the vector, \$\\mathbf{r}\$, and distance, \$|r|\$, from the node to the evaluation point

**Arguments:**
 - `node::Array{Float}` : [x y] position of node
 - `point::Array{Float}` : [x y] position of point.

**Returns**
 - `r::Vector{Float}` : vector from node to evaluation point
 - `rmag::Float` : length of panel between node and evaluation point
"""
function get_r(node, point)

    # Need to make adjustments for sqrt(0) cases
    if isapprox(point, node)
        TF = eltype(node)
        r = zeros(TF, 2)
        rmag = TF(0.0)

        return r, rmag

    else
        # Calculate vector
        r = point .- node

        # Calculate magnitude
        rmag = sqrt(r[1]^2 + r[2]^2)

        return r, rmag
    end
end

"""
    get_d(node1, node2)

Calculate panel length (between adjacent nodes).


**Arguments:**
 - `node1::Array{Float}(2)` : [x y] location of first node
 - `node2::Array{Float}(2)` : [x y] location of second node

**Returns**
 - `d::Vector{Float}` : vector from node1 to node2
 - `dmag::Float` : length of panel between node1 and node2
"""
function get_d(edges)

    # simply call get_r, since it`s exactly what is needed
    return get_r(edges[1, :], edges[2, :])
end

"""
    get_d(node1, node2)
Calculate panel length (between adjacent nodes).
**Arguments:**
 - `node1::Array{Float}(2)` : [x y] location of first node
 - `node2::Array{Float}(2)` : [x y] location of second node
**Returns**
 - `d::Vector{Float}` : vector from node1 to node2
 - `dmag::Float` : length of panel between node1 and node2
"""
function get_d(node1, node2)

    # simply call get_r, since it`s exactly what is needed
    return get_r(node1, node2)
end

"""
    get_theta(h, a)

Get angle (in radians) between panel and vector from node1 to evaluation point.

**Arguments:**
 - `h::Float` : Distance, normal to panel, between panel and evaluation point.
 - `a::Float` : Distance, tangent to panel, between node1 and evaluation point.

"""
function get_theta(h, a)
    return atan(h, a)
end

"""
    get_theta2(h, a, dmag)

Get angle (in radians) between panel and vector from node2 to evaluation point.

**Arguments:**
 - `h::Float` : Distance, normal to panel, between panel and evaluation point.
 - `a::Float` : Distance, tangent to panel, between node1 and evaluation point.
 - `dmag::Float` : Panel lentgh.

"""
function get_theta2(h, a, dmag)
    return atan(h, a - dmag)
end

"""
    get_r_normal(r1, d, dmag)

Calculate distance from panel to evalulation point in the panel normal direction.

**Arguments:**
 - `r1::Vector{Float}` : vector from node1 to evalulation point.
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_r_normal(r1, d, dmag)

    # get unit normal to panel
    nhat = get_panel_normal(d, dmag)

    # calculate h (dot product of unit normal and r1 vector
    h = r1[1] * nhat[1] + r1[2] * nhat[2]

    return h
end

"""
    get_r_tangent(r1, d, dmag)

Calculate distance from panel to evalulation point in the panel tangent direction.

**Arguments:**
 - `r1::Vector{Float}` : vector from node1 to evalulation point.
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_r_tangent(r1, d, dmag)

    # Get unit tangent vector
    that = get_panel_tangent(d, dmag)

    # calculate a (dot product of unit tangent and r1 vector)
    a = r1[1] * that[1] + r1[2] * that[2]
    return a
end

"""
    get_panel_tangent(d, dmag)

Get unit tangent to panel.

**Arguments:**
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_panel_tangent(d, dmag)
    return (dmag == 0.0) ? [0.0; 0.0] : (d / dmag)
end

"""
    get_panel_normal(d, dmag)

Get unit normal to panel.

**Arguments:**
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_panel_normal(d, dmag)

    # get unit tangent
    that = get_panel_tangent(d, dmag)

    # use fancy trick to rotate to be unit normal
    nhat = [-that[2]; that[1]]

    return nhat
end

######################################################################
#                                                                    #
#                     AXISYMMETRIC FUNCTIONS                         #
#                                                                    #
######################################################################

"""
    get_ring_geometry(paneli, panelj)

Obtain relevant geometry associated with ring singularity influence calculations.

**Arguments:**
- `paneli::FLOWFoil.AxiSymPanel` : the ith panel (the panel being influenced).
- `panelj::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).

**Returns:**
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

**Arguments:**
- `panel::FLOWFoil.AxiSymPanel` : the jth panel (the panel doing the influencing).
- `field_point::Array{Float}` : [x;r] coordinates of the field point in question.

**Returns:**
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
