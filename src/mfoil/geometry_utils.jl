"""
    get_trailing_edge_info(nodes)

Calculate various items needed for trailing edge treatment.

**Arguments:**
 - `nodes::Array{Float,2}` : Array of [x y] locations for the airfoil nodes.

**Returns:**
 - `tdp::Float` : dot product of TE bisection and TE gap unit vectors
 - `txp::Float` : "cross product" of TE bisection and TE gap unit vectors
 - `trailing_edge_gap::Float` : TE gap distance
 - `gap_edges`
 - `dte`
 - `dtemag`
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

    # gap edges
    gap_edges = [panel_edges[end, 2, :]'; panel_edges[1, 1, :]']

    # get panel vector
    dte, dtemag = get_d(panel_edges[1, 1, :], panel_edges[end, 2, :])

    if dtemag == 0.0 || bisector == zeros(2)
        return 1.0, 0.0, 0.0, gap_edges, dte, dtemag
    else
        # - Calculate tdp - #
        # normalize to get the unit vector
        ttehat = bisector / sqrt(bisector[1]^2 + bisector[2]^2)

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
