"""
    surfacenormal(xloc::Real,x::Array{<:Real,1},y::Array{<:Real,1})
Finds x,y location and angle corresponding to surface of airfoil at specified
x-location
"""
function surfacenormal(xloc::Real,x::Array{<:Real,1},y::Array{<:Real,1})
  idx = indmin(abs(x-xloc))
  xloc = x[idx]
  yloc = y[idx]
  theta = atan2((y[idx+1]-y[idx-1]),((x[idx+1]-x[idx-1])))
  return xloc,yloc,theta
end

