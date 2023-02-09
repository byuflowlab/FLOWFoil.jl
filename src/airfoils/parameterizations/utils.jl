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
    _, N = findmin(x)

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

"""
    repanel_airfoil(x,y;N)
    repanel_airfoil(xy;N)
Takes x and y coordinates of an airfoil  and uses a cosine spaced akima spline to fill in the gaps
**Arguments**
- `x::Vector{Float64}` : vector containing the x coordinates of the airfoil
- `y::Vector{Float64}` : vector containing the y components of the airfoil
- `xy::Array{Float64,2}` : Array of x and y coordinates with X in the first column and y in the 2nd
**Keyword Arguements**
- `N::Int` : Number of data points to be returned after repaneling. Will only return odd numbers, if N is even, N+1 points will be returned.
**Returns**
- `xreturn::Vector{Float64}` : Repaneled, cosine spaced x corrdinates of the airfoil
- `yreturn::Vector{Float64}` : y coordinates of the repaneled airfoil obtained using an akima spline
- `xyreturn::Array{Float64}` : If the coordinates were input as an array, this will be returned with x in the 1st column and y in the 2nd.
"""
function repanel_airfoil(x, y; N=160)
    @assert length(x) == length(y) "X and Y vectors must be the same length"

    #First normalize the airfoil to between 0 and 1
    normalize_airfoil!(x, y)

    #let's figure out the cosine spacing.
    npoints = ceil(Int, N / 2)
    akimax = cosine_spacing(npoints)

    #now we split the top and bottom of the airfoil
    x1, x2, y1, y2 = split_upper_lower(x, y)

    #Now check and see which x and y need to be reversed
    #x has to be ascending (0-->1)

    if x1[1] > x1[end]
        x1 = reverse(x1)
        y1 = reverse(y1)
    end

    if x2[1] > x2[end]
        x2 = reverse(x2)
        y2 = reverse(y2)
    end

    #do the akima spline
    akimay1 = FLOWMath.akima(x1, y1, akimax)
    akimay2 = FLOWMath.akima(x2, y2, akimax)

    #figure out which spline is on top
    if maximum(akimay1) > maximum(akimay2) #then akimay1 is on top so I need to reverse akimay2
        yreturn = [reverse(akimay2); akimay1[2:end]]

    else #otherwise akimay2 is on top so I need to reverse akimay1
        yreturn = [reverse(akimay1); akimay2[2:end]]
    end

    xreturn = [reverse(akimax); akimax[2:end]]

    return xreturn, yreturn
end

function repanel_airfoil(coordinates; N=160)
    x = coordinates[:, 1]
    y = coordinates[:, 2]

    xpane, ypane = repanel_airfoil(x, y; N=N)

    return [xpane ypane]
end
