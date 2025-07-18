"""
    flip!(x)

Flips one dimension of airfoil coordinates.

Moves airfoil left (x) or down (y) by maximum x or y coordinate then flips about the y or x axis, respectively.

# Arguments
- `x::AbstractArray{Float}` : vector of x or y coordinates
"""
function flip!(x)
    maxx = maximum(x)
    x .-= maxx
    x .*= -1.0
    return x
end

"""
    zero_y_te!(coordinates)

Places trailing edge on the x-axis.

# Arguments
- `coordinates::Array{Float}` : Array of [x y] coordinates to be updated in place.
"""
function zero_y_te!(coordinates)
    _, idxmaxx = findmax(coordinates[:, 1])
    coordinates[:, 2] .-= coordinates[idxmaxx, 2]
    return coordinates
end

"""
    rotate_coordinates!(coordinates, angle; rotation_point=[0.0; 0.0])

Rotate coordinates clockwise about `rotation_point` by `angle` in degrees.

# Arguments
- `coordinates::Array{Float}` : Array of [x y] coordinates to be updated in place.
- `angle::Float=0.0` : Angles, in degrees, by which to rotate the coordinates clockwise (positive angle will pitch airfoil up).

# Keyword Arguments
- `rotation_point::AbstractArray{Float}=[0.0; 0.0]` : Array of [x y] position of point about which to perform rotation.
"""
function rotate_coordinates!(coordinates, angle; rotation_point=[0.0; 0.0])
    # get rotation matrix
    R = [cosd(-angle) -sind(-angle); sind(-angle) cosd(-angle)]

    # rotate and translate
    for c in eachrow(coordinates)
        c[:] .-= rotation_point
        c[:] .= R * c
        c[:] .+= rotation_point
    end

    return coordinates
end

"""
    normalize_coordinates!(coordinates)

Normalize airfoil to unit chord and shift leading edge to zero. Adjusts coordinates in place.

# Arguments
- `coordinates::AbstractArray{Float}` : Array of [x y] coordinates
"""
function normalize_coordinates!(coordinates)
    x = @view(coordinates[:, 1])
    y = @view(coordinates[:, 2])

    # get current chord length
    chord = maximum(x) - minimum(x)

    # shift to zero
    x[:] .-= minimum(x)

    # normalize chord
    x[:] ./= chord

    # scale y coordinates to match
    y[:] ./= chord

    return coordinates
end

"""
    normalize_coordinates!(x, y)

Normalize airfoil to unit chord and shift leading edge to zero. Adjusts coordinates in place.

# Arguments
- `x::Array{Float}` : Array of x coordinates
- `y::Array{Float}` : Array of y coordinates
"""
function normalize_coordinates!(x, y)

    # get current chord length
    chord = maximum(x) - minimum(x)

    # shift to zero
    x[:] .-= minimum(x)

    # normalize chord
    x[:] ./= chord

    # scale y coordinates to match
    y[:] ./= chord

    return x, y
end

"""
    position_coordinates!(coordinates, scale, angle, location)

Scale, Rotate, and Transform (in that order) airfoil coordinates.

# Arguments
- `coordinates::Array{Float}` : Array of [x y] coordinates to be updated in place.

# Keyword Arguments
- `scale::Float=1.0` : Value by which to scale coordinates.
- `angle::Float=0.0` : Angles, in degrees, by which to rotate the coordinates clockwise (positive angle will pitch airfoil up).
- `location::AbstractArray{Float}=[0.0; 0.0]` : Array of [x y] position of leading edge location.
- `rotation_point::AbstractArray{Float}=[0.0; 0.0]` : Array of [x y] position of point about which to perform rotation.
- `flipped::Bool` : flag whether to flip airfoil upside down.

# Returns
- `x::Array{Float}` : array of x-coordinates
- `y::Array{Float}` : array of y-coordinates
"""
function position_coordinates!(
    coordinates;
    scale=1.0,
    angle=0.0,
    location=[0.0; 0.0],
    rotation_point=[0.0; 0.0],
    flipped=false,
)
    #flip if needed
    if flipped
        coordinates[:, 2] .*= -1.0
    end

    # scale
    coordinates .*= scale

    # get rotation matrix
    R = [cosd(-angle) -sind(-angle); sind(-angle) cosd(-angle)]

    # rotate and translate
    for c in eachrow(coordinates)
        c[:] .-= rotation_point
        c[:] = R * c
        c[:] .+= rotation_point
        c[:] .+= location
    end

    return coordinates
end

"""
    position_coordinates!(
        coordinates::Vector{AbstractArray{TF}};
        scales=[1.0],
        angles=[0.0],
        locations=[[0.0; 0.0],],
        rotation_points=[[0.0; 0.0],],
        flipped=[false],
    ) where {TF}

Multi-airfoil version of position_coordinates!.

If keyword arguments are give as single valued vectors, the same values are used for all coordinate sets.
If keyword arguments are provided as vectors of length greater than 1, they must have the same length as the set of coordinates.
For example, if scaling 3 airfoils, there will be a vector of 3 airfoil coordinate sets input and `scales` must either be a one element vector or a vector of length 3.
"""
function position_coordinates!(
    coordinates::Vector{AbstractArray{TF}};
    scales=[1.0],
    angles=[0.0],
    locations=[[0.0; 0.0]],
    rotation_points=[[0.0; 0.0]],
    flipped=[false],
) where {TF}

    # check that things are the right lengths

    if length(scales) > 1
        @assert length(scales) == length(coordinates) "scales must either be length 1 or of the same length as coordinates"
    end

    if length(angles) > 1
        @assert length(angles) == length(coordinates) "angles must either be length 1 or of the same length as coordinates"
    end

    if length(locations) > 1
        @assert length(locations) == length(coordinates) "locations must either be length 1 or of the same length as coordinates"
    end

    if length(rotation_points) > 1
        @assert length(rotation_points) == length(coordinates) "rotation_points must either be length 1 or of the same length as coordinates"
    end

    if length(flipped) > 1
        @assert length(flipped) == length(coordinates) "flipped must either be length 1 or of the same length as coordinates"
    end

    # Loop through coordinates
    for (i, c) in enumerate(coordinates)
        position_coordinates!(
            @view(c[:, :]);
            scale=length(scales) > 1 ? scales[i] : scales,
            angle=length(angles) > 1 ? angles[i] : angles,
            location=length(locations) > 1 ? locations[i] : locations,
            rotation_point=if length(rotation_points) > 1
                rotation_points[i]
            else
                rotation_points
            end,
            flipped=length(flipped) > 1 ? flipped[i] : flipped,
        )
    end

    return coordinates
end
