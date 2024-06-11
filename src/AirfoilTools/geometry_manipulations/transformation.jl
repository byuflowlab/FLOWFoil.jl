"""
    flipx(x::Array{<:Real,1})
Flips airfoil x coordinates about y axis
"""
function flipx(x::Array{<:Real,1},y::Array{<:Real,1})
  maxx = maximum(x)
  x-=maxx
  x=-x
  return x,y
end

"""
    zerochordy(x::Array{<:Real,1},y::Array{<:Real,1})
Zeros chord line after rotation
"""
function zerochordy(x::Array{<:Real,1},y::Array{<:Real,1})
  idxmaxx = indmax(x)
  y-=y[idxmaxx]
  return x,y
end

"""
    rotateaf(x::Array{<:Real,1},y::Array{<:Real,1},angle)
Rotates airfoil ccw about origin
"""
function rotateaf(x::Array{<:Real,1},y::Array{<:Real,1},angle)
  x = cos(angle)*x-sin(angle)*y
  y = cos(angle)*y+sin(angle)*x
  return x,y
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
    transformaf(x::Array{<:Real,1},y::Array{<:Real,1},scale::Real,
    twist::Real,loc::Array{<:Real,1})
Scales, rotates, and translates an airfoil.
"""
function transformaf(x::Array{<:Real,1},y::Array{<:Real,1},scale::Real,twist::Real,loc::Array{<:Real,1})
  # Rotation matrix
  R = [[cos(twist) 0 sin(twist)]
      [0 1 0]
      [-sin(twist) 0 cos(twist)]]

  airfoil = zeros(Real,length(x),3)
  for i = 1:length(x)
    airfoil[i,1] = x*scale
    airfoil[i,2] = 0.0
    airfoil[i,3] = y*scale
    airfoil[i,:] = R*airfoil[i,:] + loc
  end
  return airfoil
end

"""
    transformaf(afdata::Array{<:Real,3},scale::Array{<:Real,1},
    twist::Array{<:Real,1},loc::Array{<:Real,2})
Multi-airfoil implementation of transformairfoil. First dimension of each
input/output corresponds to specific airfoil.
"""
function transformaf(afdata::Array{<:Real,3},scale::Array{<:Real,1},twist::Array{<:Real,1},loc::Array{<:Real,2})
  airfoils = zeros(Real,size(afdata,1),size(afdata,2),3)
  for i = 1:size(afdata,1)
    airfoil = transformairfoil(afdata[i,:,1],afdata[i,:,2],scale[i],twist[i],loc[i,:])
    for j = 1:size(afdata,2)
      for k = 1:size(afdata,3)
        airfoils[i,j,k] = airfoil[j,k]
      end
    end
  end
  return airfoils
end


