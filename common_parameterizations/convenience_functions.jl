#=
Various convenience functions for dealing with airfoil parameterizations and coordinates.

=#

"""
    function cosinespacing(N=180)

Calculate N cosine spaced points.

**Arguments:**
 - 'N::Int' : Number of points

"""
function cosinespacing(N=80)
    return [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]
end

"""
    split_upper_lower(x, z)

Split the upper and lower halves of the airfoil coordinates. Assumes odd number of coordinates (leading edge repeated).

**Arguments:**
 - 'x::Array{Float}' : Array of x coordinates
 - 'z::Array{Float}' : Array of z coordinates

**Returns:**
 - 'xu::Array{Float}' : Array of upper half of x coordinates
 - 'xl::Array{Float}' : Array of lower half of x coordinates
 - 'zu::Array{Float}' : Array of upper half of z coordinates
 - 'zl::Array{Float}' : Array of lower half of z coordinates

"""
function split_upper_lower(x, z)

    # get half length of geometry coordinates
    N = Int(ceil(length(x) / 2.0))

    return x[1:N], x[N:end], z[1:N], z[N:end]
end

"""
    normalize_airfoil!(x, z)

Normalize airfoil to unit chord and shift leading edge to zero. Adjusts coordinates in place.

**Arguments:**
 - 'x::Array{Float}' : Array of x coordinates
 - 'z::Array{Float}' : Array of z coordinates

"""
function normalize_airfoil!(x, z)
    chord = maximum(x) - minimum(x) #get current chord length
    x .-= minimum(x) #shift to zero
    x ./= chord #normalize chord
    z ./= chord #scale z coordinates to match

    return nothing
end
