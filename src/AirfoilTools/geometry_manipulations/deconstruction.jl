"""
    split_upper_lower(x, z; idx::Integer=nothing)

Split the upper and lower halves of the airfoil coordinates.

Assumes leading edge point is at first minimum x value if `idx` is not provided.
Returns the upper and lower coordinates each with the leading edge point.

# Arguments:
 - `x::Vector{Float}` : Vector of x coordinates
 - `z::Vector{Float}` : Vector of z coordinates

# Keyword Arguments:
 - `idx::Integer` : optional index at which to split the coordinates

# Returns:
 - `xu::Vector{Float}` : Vector of upper half of x coordinates
 - `xl::Vector{Float}` : Vector of lower half of x coordinates
 - `zu::Vector{Float}` : Vector of upper half of z coordinates
 - `zl::Vector{Float}` : Vector of lower half of z coordinates

"""
function split_upper_lower(x, z; idx=nothing)

    # get half length of geometry coordinates
    if isnothing(idx)
        _, idx = findmin(x)
    end

    return x[1:idx], x[idx:end], z[1:idx], z[idx:end]
end
