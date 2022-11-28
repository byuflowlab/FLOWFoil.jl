#=
Collection of various airfoil geometry manipulation functions.

Authors: Taylor McDonnell

=#

module AirfoilManip

export flipx,zerochordy,rotateaf,splitaf,transformaf,addtepoints,fliplevelfix,
  fixaf,naca4,surfacenormal,cosinespace

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

"""
    addtepoints(x::Array{<:Real,1},y::Array{<:Real,1})
Adds points along the trailing edge of an airfoil
"""
function addtepoints(x::Array{<:Real,1},y::Array{<:Real,1})
  # Find appropriate trailing edge spacing
  # Take average of spacing of points adjacent to trailing edge.
  teapproxspacing = (sqrt((x[end-1]-x[end])^2+(y[end-1]-y[end])^2)+
    sqrt((x[2]-x[1])^2+(y[2]-y[1])^2))/2
  # Find size of trailing edge
  tesize = abs(y[end]-y[1])
  # See how many points is appropriate along trailing edge
  N = round(Integer,tesize/teapproxspacing)
  x,y = addtepoints(x,y,N)
  return x,y
end

"""
    addtepoints(x::Array{<:Real,1},y::Array{<:Real,1},N::Integer)
Adds points along the trailing edge of an airfoil
"""
function addtepoints(x::Array{<:Real,1},y::Array{<:Real,1},N::Integer)
  x = vcat(x[1:(end-1)],ones(N-1))
  yte = linspace(y[end],y[1],N)
  y = vcat(y[1:(end-1)],yte[1:(end-1)])
  return x,y
end

"""
    fliplevelfix(x::Array{<:Real,1},y::Array{<:Real,1},angle,plot=true)
Combines fixaf() rotate() fixaf() zerochordy() flipx() to level, flip,
and normalize an airfoil for a given angle
"""
function fliplevelfix(x::Array{<:Real,1},y::Array{<:Real,1},angle::Real)
  x,y = fixaf(x,y)
  x,y = rotateaf(x,y,angle)
  x,y = fixaf(x,y)
  x,y = flipx(x,y)
  x,y = zerochordy(x,y)
  return x,y
end

"""
    fixaf(x,y)
Adjusts coordinates of airfoil to loop from bottom edge trailing edge to top
trailing edge
"""
function fixaf(x::Array{<:Real,1},y::Array{<:Real,1})
  # Ensure proper scaling
  minx = minimum(x)
  x = x-minx
  maxx = maximum(x)
  x = x/maxx
  y = y/maxx

  # Check if only top or bottom specified.  Assume symmetric and reflect if true.
  if ((x[1] == 0.0) && (x[end] == 1.0)) ||
    ((x[1] == 1.0) && (x[end] == 0.0))
    idx = find(y.!=0.0)
    xsym = x[idx]
    ysym = -y[idx]
    xsym = flipdim(xsym,1)
    ysym = flipdim(ysym,1)
    x = vcat(x,xsym)
    y = vcat(y,ysym)
  end

  # Find trailing edge points
  idxte = indmax(x)

  #blunt trailing edge if trailing edge is at angle
  if length(idxte) == 1 && y[idxte] !=0.0
    blunt = true
    # check which direction is more vertical
    coord = hcat(x,y)
    if idxte == length(x)
      dotp1idx=1
    else
      dotp1idx=idxte+1
    end
    dotp1 = abs(dot(coord[idxte,:]-coord[dotp1idx,:],[0.0,1.0]))/(norm(coord[idxte,:]-coord[dotp1idx,:])*norm([0.0,1.0]))
    if idxte == 1
      dotm1idx=length(x)
    else
      dotm1idx=idxte-1
    end
    dotm1 = abs(dot(coord[idxte,:]-coord[dotm1idx,:],[0.0,1.0]))/(norm(coord[idxte,:]-coord[dotm1idx,:])*norm([0.0,1.0]))
    if dotp1 > dotm1
      idxte = [idxte,dotp1idx]
    else
      idxte = [idxte,dotm1idx]
    end
  elseif length(idxte) > 1
    blunt = true
  else
    blunt = false
  end

  idxstart = idxte[indmin(y[idxte])]
  # If trailing edge of pressure side is not first element, change accordingly
  if (idxstart !=1)
    x = vcat(x[idxstart:end],x[1:(idxstart-1)])
    y = vcat(y[idxstart:end],y[1:(idxstart-1)])
    idxte = 1
  end
  # Find leading edge
  idxle = indmin(x)
  # Check if not clockwise loop, flip if necessary
  if (sum(y[1:idxle]) > sum(y[idxle:end]))
    x = flipdim(x,1)
    y = flipdim(y,1)
    idxte = length(x)
  end
  # Check if blunt trailing edge and repeat te if not
  if !blunt
    x = vcat(x,x[idxte])
    y = vcat(y,y[idxte])
  end
  return x,y
end

"""
    naca4(naca::String,N::Integer=200)
Constructs a standard NACA 4 digit airfoil with N cosine spaced points
"""
function naca4(naca::String,N::Integer=200)
    return naca4(naca,cosinespace(N))
end

"""
    naca4(naca::String,xc::Array{<:Real,1})
Constructs a standard NACA 4 digit airfoil with normalized x coordinates
specified in xc.
"""
function naca4(naca::String,xc::Array{<:Real,1})
  # Parse identifier
  m = parse(naca[1:1])/100.0
  p = parse(naca[2:2])/10.0
  tt = parse(naca[3:4])/100.0

  # Find camber line if not symmetric
  yc = zeros(Real,length(x))
  theta = zeros(Real,length(x))
  if (m != 0.0) && (p != 0.0)
    for i = 1:length(x)
      if (x[i] <= p)
        yc[i] = m/p^2*(2*p*x[i]-(x[i])^2.0)
        theta[i] = atan(2*m/p^2*(p-x[i]))
      else
        yc[i] = m/(1-p)^2*((1-2*p)+2*p*x[i]-(x[i])^2.0)
        theta[i] = atan(2*m/(1-p)^2*(p-x[i]))
      end
    end
  end

  # Get thickness
  yt = zeros(Real,length(x))
  for i = 1:length(x)
    yt[i] = tt/0.2*(0.2969*sqrt(x[i])-0.1260*x[i]-0.3516*x[i]^2+0.2843*x[i]^3-0.1015*x[i]^4)
  end

  # Combine thickness and camber
  leidx = indmin(x)
  y = zeros(Real,length(x))
  for i = 1:length(x)
    if i < leidx
      x[i] = x[i]+yt[i]*sin(theta[i])
      y[i] = yc[i]-yt[i]*cos(theta[i])
    else
      x[i] = x[i]-yt[i]*sin(theta[i])
      y[i] = yc[i]+yt[i]*cos(theta[i])
    end
  end
  x,y = fixaf(x,y)
  return x,y
end

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

"""
    cosinespace(N::Integer)
Returns cosine spaced x coordinates
"""
function cosinespace(N::Integer)
    x = zeros(Real,N+1)
    for i=1:N+1
      zeta = 2*pi/(N)*(i-1)
      x[i] = 0.5*(cos(zeta)+1)
    end
    return x
end

end # module
