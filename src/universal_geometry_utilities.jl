"""
    get_d(node1, node2)

Calculate panel length (between adjacent nodes).


# Arguments:
 - `node1::Array{Float}(2)` : [x y] location of first node
 - `node2::Array{Float}(2)` : [x y] location of second node

# Returns
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
# Arguments:
 - `node1::Array{Float}(2)` : [x y] location of first node
 - `node2::Array{Float}(2)` : [x y] location of second node
# Returns
 - `d::Vector{Float}` : vector from node1 to node2
 - `dmag::Float` : length of panel between node1 and node2
"""
function get_d(node1, node2)

    # simply call get_r, since it`s exactly what is needed
    return get_r(node1, node2)
end

"""
    function get_r(node,point)

Calculate the vector, \$\\mathbf{r}\$, and distance, \$|r|\$, from the node to the evaluation point

# Arguments:
 - `node::Array{Float}` : [x y] position of node
 - `point::Array{Float}` : [x y] position of point.

# Returns
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
    get_panel_tangent(d, dmag)

Get unit tangent to panel.

# Arguments:
 - `d::Vector{Float}` : vector from node1 to node2.
 - `dmag::Float` : panel length

"""
function get_panel_tangent(d, dmag)
    return (dmag == 0.0) ? [0.0; 0.0] : (d / dmag)
end

"""
    get_panel_normal(d, dmag)

Get unit normal to panel.

# Arguments:
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
