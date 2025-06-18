# Airfoil Geometry Manipulation Tools

Here we include the variety of methods for manipulating airfoil geometries in useful ways implemented in AirfoilTools.

## Deconstruction

It is often convenient to deconstruct an airfoil into its upper and lower halves.  The `split_upper_lower` function makes this process straightforward.

```@example split
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

x, y = naca4()

xl, xu, yl, yu = split_upper_lower(x, y)

plot(xl, yl; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Lower Side")
plot!(xu, yu; label="Upper Side")
```

```@docs
AirfoilTools.split_upper_lower
```

## Translation

There are various transformation functions that can also be helpful in various situations.

We can flip an airfoil.
```@example trans
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

x, y = naca4()

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Nominal")
plot!(flip!(copy(x)), y; label="x Flipped")
```

!!! note
    Note that this function both flips and translates. It's primarily useful for flipping the x coordinates.  If you want to flip the y coordinate, applying a simple negative would be best

```@example trans
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Nominal")
plot!(x, flip!(copy(y)); label="y Flipped")
plot!(x, -y; label="y Negated")
```

Another thing we could do is translate the airfoil so that the trailing edge is at y=0.

```@example trans
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

x, y = naca4()
y .+= 0.2

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Nominal")

xy = [x y]
zero_y_te!(xy)
plot!(xy[:,1], xy[:,2]; label="TE zeroed")
```

!!! note
    This only vertically translates the geometry, it does not rotate things to put the leading edge on the axis as well if the geometry is rotated.


```@docs
AirfoilTools.zero_y_te!
```

## Rotation

We can also rotate airfoils about arbitrary points.

```@example rotate
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

x, y = naca4()

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Nominal")

xy = [x y]
angle = 10.0 # degrees

# rotate by angle about default point: (0,0)
rotate_coordinates!(xy, angle)

plot!(xy[:,1], xy[:,2], label="rotated about (0,0)")

# rotate again but about the trailing edge
xy2 = [x y]
rotate_coordinates!(xy2, angle; rotation_point=[1.0,0.0])

plot!(xy2[:,1], xy2[:,2], label="rotated about (1,0)")
```

```@docs
AirfoilTools.rotate_coordinates!
```

## Normalization

It can also be convenient to normalize coordinates to have a chord length of one.

```@example norm
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

x, y = 2.0.*naca4()

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Nominal")

normalize_coordinates!(x, y)

plot!(x, y; label="Normalized")
```

!!! note
    This function is designed to go from a nominal airfoil (length 1, leading edge and trailing edges on the axis, etc.) to something else.  The operations are in the order: scale, rotate, then translate.

```@docs
AirfoilTools.normalize_coordinates!
```

## Scale, Rotate, and Translate

Sometimes, we might want do several operations together, which we can with the `position_coordinates!` function.

```@example pose
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

x, y = naca4()
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Unscaled")

xy = [x y]
position_coordinates!(
    xy; scale=0.8, angle=3.0, location=[-0.2; 0.1], rotation_point=[0.0; 0.0], flipped=true
)

plot!(xy[:, 1], xy[:, 2]; label="Re-positioned")
```

```@docs
AirfoilTools.position_coordinates!
```

## Re-definition

Sometimes, we may want to "re-panel" an airfoil, (a good example is the PANE command in Xfoil).  The `repanel_airfoil` method splines the provided geometry and then resamples it using cosine spacing to give a higher density of coordinates at the leading and trailing edges.

```@example repanel
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

x, y = naca4(; x = range(0.0,1.0,30))
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Linear Spaced", marker=true, markersize=4)

x_rp, y_rp = repanel_airfoil(x,y; N=161)

scatter!(x_rp, y_rp; label="Cosine Spaced", markersize=1.5)
```

!!! note
    Note that if you simply have too few coordinates to begin with, repaneling the airfoil isn't going to smooth out an inaccurate leading edge.  In such a case you should consider fitting the geometry with one of the airfoil parameterizazation methods, which will ensure a round leading edge.

```@docs
AirfoilTools.repanel_airfoil
AirfoilTools.repanel_revolution
```

## Contributing

We welcome the addition of more convenience functions for airfoil geometry manipulation.
