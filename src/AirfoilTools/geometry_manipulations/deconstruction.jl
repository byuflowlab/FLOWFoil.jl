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
    splitaf(x,y)
Splits airfoil into x and y coordinates of pressure side and suction side
"""
function splitaf(x::Array{<:Real,1},y::Array{<:Real,1})
  lei = indmin(abs(x))
  xps,xss,yps,yss = splitairfoil(x,y,lei)
  return xps,xss,yps,yss
end

"""
    splitaf(x,y,idx)
Splits airfoil into x and y coordinates at index location idx
"""
function splitaf(x::Array{<:Real,1},y::Array{<:Real,1},idx::Integer)
  x1 = x[1:idx]
  y1 = y[1:idx]
  x2 = x[idx:end]
  y2 = y[idx:end]
  return x1,x2,y1,y2
end


