

"""
"""
function get_vtbar(theta1, theta2)
    return (theta2 - theta1) / (2 * pi)
end

"""
"""
function get_vttilde(theta1, theta2, ln1, ln2, dmag, h, a)
    return (h * (ln2 - ln1) + a * (theta2 - theta1)) / (2 * pi * dmag)
end

"""
"""
function get_vnbar(ln1, ln2)
    return (ln2 - ln1) / (2 * pi)
end

"""
"""
function get_vntilde(theta1, theta2, ln1, ln2, dmag, h, a)
    return (a * (ln2 - ln1) + d - h * (theta2 - theta1)) / (2 * pi * d)
end

"""
"""
function vt(mesh, gamma)

    # get system size
    N, Ns = FLOWFoil.sizesystem(mesh)

    # get nodes for convenience
    nodes = mesh.meshes[1].airfoil_nodes

    #initialize NxN coefficient matrix
    vt = Vector{Float64}(undef, N)

    vtbar = get_vtbar(theta1, theta2)
    return vttilde = get_vttilde(theta1, theta2, ln1, ln2, dmag, h, a)
end

