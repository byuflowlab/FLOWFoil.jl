"""
    split_upper_lower(x, z; idx::Integer=nothing)

Split the upper and lower halves of the airfoil coordinates.

Assumes leading edge point is at first minimum x value if `idx` is not provided.
Returns the upper and lower coordinates each with the leading edge point.
Assumes airfoil is defined clockwise starting at the trailing edge.

# Arguments:
 - `x::AbstractArray{Float}` : Vector of x coordinates
 - `z::AbstractArray{Float}` : Vector of z coordinates

# Keyword Arguments:
 - `idx::Integer` : optional index at which to split the coordinates

# Returns:
 - `xl::AbstractArray{Float}` : Vector of lower half of x coordinates
 - `xu::AbstractArray{Float}` : Vector of upper half of x coordinates
 - `zl::AbstractArray{Float}` : Vector of lower half of z coordinates
 - `zu::AbstractArray{Float}` : Vector of upper half of z coordinates

"""
function split_upper_lower(x, z; idx=nothing)

    # get half length of geometry coordinates
    if isnothing(idx)
        _, idx = findmin(x)
    end

    return x[1:idx], x[idx:end], z[1:idx], z[idx:end]
end

"""
    split_upper_lower(coordaintes; idx::Integer=nothing)

Split the upper and lower halves of the airfoil coordinates.

Assumes leading edge point is at first minimum x value if `idx` is not provided.
Returns the upper and lower coordinates each with the leading edge point.
Assumes airfoil is defined clockwise starting at the trailing edge.

# Arguments:
 - `coordinates::Matrix{Float}` : Matrix of [x z] coordinates

# Keyword Arguments:
 - `idx::Integer` : optional index at which to split the coordinates

# Returns:
 - `xl::AbstractArray{Float}` : View of lower half of x coordinates
 - `xu::AbstractArray{Float}` : View of upper half of x coordinates
 - `zl::AbstractArray{Float}` : View of lower half of z coordinates
 - `zu::AbstractArray{Float}` : View of upper half of z coordinates

"""
function split_upper_lower(coordinates; idx=nothing)
    x = @view(coordinates[:, 1])
    z = @view(coordinates[:, 2])

    # get half length of geometry coordinates
    if isnothing(idx)
        _, idx = findmin(x)
    end

    return x[1:idx], x[idx:end], z[1:idx], z[idx:end]
end
