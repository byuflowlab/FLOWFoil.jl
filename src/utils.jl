
"""
    linear_transform(range1, range2, values)

Linear transfrom of values from range (source_range[1], raend) to (target_range[1], target_range[end])

**Arguments:**
- `source_range::Vector{Float{` : range values come from
- `target_range::Vector{Float}` : range onto which we are transforming
- `source_values::Array{Float}` : array of source_values to transform

**Returns:**
 - `target_values::Array{Float}` : array of transformed source_values onto target range
"""
function linear_transform(source_range, target_range, source_values)
    return target_range[1] .+
           (target_range[end] - target_range[1]) .* (source_values .- source_range[1]) /
           (source_range[end] - source_range[1])
end

"""
    smooth_distributions(::Order, surface_location, surface_values, npanels) end

Generates smooth surface distribution values.

**Arguments:**
- `o::Order` : Order of input geometry
- `surface_location::Array{Float}` : location where distribution values lie.
- `surface_values::Array{Float}` : surface distribution values
- `npanels::Int` : Number of panels to use on the top and bottom surface for smoothing (total panels = 2*npanels-1)

**Returns:**
- `distribution::Vector{Float}` : Smoothed surface distribution
- `xsmooth::Vector{Float}` : Smoothed x-coordinates

# TODO: Consider returning spline objects rather than another set of discrete values.
"""
function smooth_distributions(::Order, surface_location, surface_values, npanels) end

function smooth_distributions(::Linear, panel_edges, distribution, npanels)

    #= NOTE:
        Akima splines in FLOWMath require the 'x' values to be monotonically ascending.
        Therefore, we need to get all the panel edge points and then divide them into top and bottom in order to create our splines.
    =#

    # - Get 'x' values from panel edges - #
    x = [panel_edges[:, 1, 1]; panel_edges[end, 2, 1]]

    # - Split the 'x' values - #

    # find the minimum and index
    minx, minidx = findmin(x)

    # the bottom needs to be flipped to ascend monotonically
    xbot = x[minidx:-1:1]
    # the top is already in the right direction
    xtop = x[minidx:end]

    # - Get smooth 'x' values from cosine spacing - #
    # Get cosine spaced values from zero to one.
    xcosine = cosine_spacing(npanels)

    # - Transform the cosine spaced values to the minimum and maximum points - #
    # Get the maximum x value
    maxx = maximum(x)

    xsmooth = linear_transform([0.0; 1.0], [minx; maxx], xcosine)

    # - Generate smooth distribution - #
    distbot = FLOWMath.akima(xbot, distribution[minidx:-1:1], xsmooth)
    disttop = FLOWMath.akima(xtop, distribution[minidx:end], xsmooth)

    # - Combine distribution and x values - #
    xs = [reverse(xsmooth); xsmooth[2:end]]
    dist = [reverse(distbot); disttop[2:end]]

    # - Return - #
    return dist, xs
end

function smooth_distributions(::Constant, panel_center, distribution, npanels)

    #= NOTE:
        Akima splines in FLOWMath require the 'x' values to be monotonically ascending.
        Therefore, we need to get all the panel edge points and then divide them into top and bottom in order to create our splines.
    =#

    # - Get 'x' values from panel centers - #
    x = panel_center[:, 1]

    # - Split the 'x' values - #

    # find the minimum and index
    minx, minidx = findmin(x)

    # the bottom needs to be flipped to ascend monotonically
    xbot = x[minidx:-1:1]
    # the top is already in the right direction
    xtop = x[minidx:end]

    # - Get smooth 'x' values from cosine spacing - #
    # Get cosine spaced values from zero to one.
    xcosine = cosine_spacing(npanels)

    # - Transform the cosine spaced values to the minimum and maximum points - #
    # Get the maximum x value
    maxx = maximum(x)

    xsmooth = linear_transform([0.0; 1.0], [minx; maxx], xcosine)

    # - Generate smooth distribution - #
    distbot = FLOWMath.akima(xbot, distribution[minidx:-1:1], xsmooth)
    disttop = FLOWMath.akima(xtop, distribution[minidx:end], xsmooth)

    # - Combine distribution and x values - #
    xs = [reverse(xsmooth); xsmooth[2:end]]
    dist = [reverse(distbot); disttop[2:end]]

    # - Return - #
    return dist, xs
end
