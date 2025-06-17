"""
    split_upper_lower(x, y; idx::Integer=nothing)

Split the upper and lower halves of the airfoil coordinates.

Assumes leading edge point is at first minimum x value if `idx` is not provided.
Returns the upper and lower coordinates each with the leading edge point.
Assumes airfoil is defined clockwise starting at the trailing edge.

# Arguments
 - `x::AbstractArray{Float}` : Vector of x coordinates
 - `y::AbstractArray{Float}` : Vector of y coordinates

# Keyword Arguments
 - `idx::Integer` : optional index at which to split the coordinates

# Returns
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates
 - `yl::AbstractArray{Float}` : Vector of lower half of y coordinates
 - `yu::AbstractArray{Float}` : Vector of upper half of y coordinates
"""
function split_upper_lower(x, y; idx=nothing)

    # get half length of geometry coordinates
    if isnothing(idx)
        _, idx = findmin(x)
    end

    return x[1:idx], x[idx:end], y[1:idx], y[idx:end]
end

"""
    split_upper_lower(coordaintes; idx::Integer=nothing)

Split the upper and lower halves of the airfoil coordinates.

Assumes leading edge point is at first minimum x value if `idx` is not provided.
Returns the upper and lower coordinates each with the leading edge point.
Assumes airfoil is defined clockwise starting at the trailing edge.

# Arguments
 - `coordinates::Matrix{Float}` : Matrix of [x y] coordinates

# Keyword Arguments
 - `idx::Integer` : optional index at which to split the coordinates

# Returns
 - `xl::AbstractArray{Float}` : View of lower half of x coordinates
 - `xu::AbstractArray{Float}` : View of upper half of x coordinates
 - `yl::AbstractArray{Float}` : View of lower half of y coordinates
 - `yu::AbstractArray{Float}` : View of upper half of y coordinates
"""
function split_upper_lower(coordinates; idx=nothing)
    x = @view(coordinates[:, 1])
    y = @view(coordinates[:, 2])

    # get half length of geometry coordinates
    if isnothing(idx)
        _, idx = findmin(x)
    end

    return x[1:idx], x[idx:end], y[1:idx], y[idx:end]
end
