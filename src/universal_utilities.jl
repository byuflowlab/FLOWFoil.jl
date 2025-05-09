#= NOTE:
    Used in problem definition function to help count number of bodies, the coordinates of which are a tuple of vectors if multiple bodies are being analyzed together.
=#
import Base.size
function size(t::Tuple)
    return length(t)
end

"""
    linear_transform(range1, range2, values)

Linear transfrom of values from range (source_range[1], raend) to (target_range[1], target_range[end])

# Arguments:
- `source_range::Vector{Float{` : range values come from
- `target_range::Vector{Float}` : range onto which we are transforming
- `source_values::Array{Float}` : array of source_values to transform

# Returns:
 - `target_values::Array{Float}` : array of transformed source_values onto target range
"""
function linear_transform(source_range, target_range, source_values)
    return target_range[1] .+
           (target_range[end] - target_range[1]) .* (source_values .- source_range[1]) /
           (source_range[end] - source_range[1])
end

"""
    smooth_distributions(method, panel_geometry, surface_values, npoints)

Generates smooth surface distribution values.

# Arguments:
- `panel_geometry::NamedTuple` : NamedTuple that comes from `generate_panel_geometry`. Must contain `panel_edges`.
- `distribution::Array{Float}` : surface distribution values
- `npoints::Int` : Total number of points in the new distribution covering the entire airfoil. (There will be `(npoints+1)/2` points on the top and bottom surfaces, respectively.)

# Returns:
- `distribution::Vector{Float}` : Smoothed surface distribution
- `xsmooth::Vector{Float}` : Smoothed x-coordinates

# TODO: Consider returning spline objects rather than another set of discrete values.
"""
function smooth_distributions(method::Mfoil, panel_geometry, distribution, npoints)

    #= NOTE:
        Akima splines in FLOWMath require the 'x' values to be monotonically ascending.
        Therefore, we need to get all the panel edge points and then divide them into top and bottom in order to create our splines.
    =#

    # get panel_edges from panel_geometry
    panel_edges = panel_geometry.panel_edges

    # - Get 'x' values from panel edges - #
    x = [panel_edges[:, 1, 1]; panel_edges[end, 2, 1]]

    # - Split the 'x' values - #

    # find the minimum and index
    minidx = findall(==(minimum(x)), x)
    idxbot = minidx[1]
    if length(minidx) > 1
        idxtop = minidx[end]
    else
        idxtop = minidx[1]
    end

    # the bottom needs to be flipped to ascend monotonically
    xbot = x[idxbot:-1:1]
    # the top is already in the right direction
    xtop = x[idxtop:end]

    # - Get smooth 'x' values from cosine spacing - #
    # Get cosine spaced values from zero to one.
    xcosine = at.split_cosine_spacing((npoints + 1) / 2)

    # - Transform the cosine spaced values to the minimum and maximum points - #
    # Get the maximum x value
    minx = minimum(x)
    maxx = maximum(x)

    xsmooth = linear_transform([0.0; 1.0], [minx; maxx], xcosine)

    # - Generate smooth distribution - #
    distbot = FLOWMath.akima(xbot, distribution[idxbot:-1:1], xsmooth)
    disttop = FLOWMath.akima(xtop, distribution[idxtop:end], xsmooth)

    # - Combine distribution and x values - #
    xs = [reverse(xsmooth); xsmooth[2:end]]
    dist = [reverse(distbot); disttop[2:end]]

    # - Return - #
    return dist, xs
end

function smooth_distributions(
    method::Lewis, panel_geometry, distribution, npoints; body_of_revolution=false
)

    #= NOTE:
        Akima splines in FLOWMath require the 'x' values to be monotonically ascending.
        Therefore, we need to get all the panel edge points and then divide them into top and bottom in order to create our splines.
    =#

    panel_center = panel_geometry.panel_center

    # - Get 'x' values from panel centers - #
    x = panel_center[:, 1]

    # - Split the 'x' values - #

    # find the minimum and index
    minidx = findall(==(minimum(x)), x)
    idxbot = minidx[1]
    if length(minidx) > 1
        idxtop = minidx[end]
    else
        idxtop = minidx[1]
    end

    # the top is already in the right direction
    xtop = x[idxtop:end]

    # - Get smooth 'x' values from cosine spacing - #
    # Get cosine spaced values from zero to one.
    if body_of_revolution
        xcosine = cosine_spacing(npoints)
    else
        xcosine = cosine_spacing((npoints + 1) / 2)
    end

    # - Transform the cosine spaced values to the minimum and maximum points - #
    # Get the maximum x value
    minx = minimum(x)
    maxx = maximum(x)

    xsmooth = linear_transform([0.0; 1.0], [minx; maxx], xcosine)

    # - Generate smooth distribution - #
    disttop = FLOWMath.akima(xtop, distribution[idxtop:end], xsmooth)

    # - Get the bottom surface stuff if not a body of revolution - #
    if !body_of_revolution
        # the bottom needs to be flipped to ascend monotonically
        xbot = x[idxbot:-1:1]
        distbot = FLOWMath.akima(xbot, distribution[idxbot:-1:1], xsmooth)
        # - Combine distribution and x values - #
        xs = [reverse(xsmooth); xsmooth[2:end]]
        dist = [reverse(distbot); disttop[2:end]]
    else
        xs = xsmooth
        dist = disttop
    end

    # - Return - #
    return dist, xs
end

"""
    dot(A, B) = sum(a * b for (a, b) in zip(A, B))

A faster dot product.
"""
dot(A, B) = sum(a * b for (a, b) in zip(A, B))
