#=
Various convenience functions for dealing with airfoil parameterizations and coordinates.
TODO: move all of these to the airfoil manipulations file.

=#

"""
    function cosine_spacing(N=180)

Calculate N cosine spaced points.

**Arguments:**
 - `N::Int` : Number of points

"""
function cosine_spacing(N=80)
    return [0.5 * (1 - cos(pi * (i - 1) / (N - 1))) for i in 1:N]
end

"""
    split_upper_lower(x, z)

Split the upper and lower halves of the airfoil coordinates. Assumes odd number of coordinates (leading edge repeated).

**Arguments:**
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates

**Returns:**
 - `xu::Array{Float}` : Array of upper half of x coordinates
 - `xl::Array{Float}` : Array of lower half of x coordinates
 - `zu::Array{Float}` : Array of upper half of z coordinates
 - `zl::Array{Float}` : Array of lower half of z coordinates

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
 - `x::Array{Float}` : Array of x coordinates
 - `z::Array{Float}` : Array of z coordinates

"""
function normalize_airfoil!(x, z)
    chord = maximum(x) - minimum(x) #get current chord length
    x .-= minimum(x) #shift to zero
    x ./= chord #normalize chord
    z ./= chord #scale z coordinates to match

    return nothing
end

"""
    position_coordinates!(meshes, scales, angles, locations)

Take in meshes and adjust scale, leading edge location, and angle of attack of the individual meshes in the system.  Updates mesh objects in place.

**Arguments:**
 - `meshes::Array{PlanarMesh}` : Array of mesh objects.
 - `scales::Array{Float}` : Array of numbers by which to scale respective meshes.
 - `angles::Array{Float}` : Array of angles, in degrees, by which to rotate respective meshes (positive = pitch up).
 - `locations::Array{Array{Float}}` : Array of [x y] positions of leading edges for respective meshes.

**Keyword Arguments:**
- `flipped::Bool` : flag whether to flip airfoil upside down

**Returns:**
- `xcoordinates::Array{Float}` : array of x-coordinates
- `zcoordinates::Array{Float}` : array of z-coordinates

"""
function position_coordinates(
    coordinates, scale, angle, location; flipped=false, constant_point=[0.0 0.0]
)

    #flip if needed
    if flipped
        coordinates[:, 2] .*= -1.0
        reverse!(coordinates; dims=1)
    end

    # scale
    coordinates .*= scale

    coordinates .-= constant_point

    # get rotation matrix
    R = [cosd(-angle) -sind(-angle); sind(-angle) cosd(-angle)]

    # rotate and translate
    for j in 1:length(coordinates[:, 1])
        coordinates[j, :] = R * coordinates[j, :]
        coordinates[j, :] .+= location
    end
    coordinates .+= constant_point

    return coordinates[:, 1], coordinates[:, 2]
end
