# Common Parameterizations

### Parameter Types
Most of the parameterizations below have associated composite types whose fields are the parameters used in defining the airfoil geometries.  Each of these composite types is defined using the `@kwdef` macro such that the user does not need to remember the order of the fields, but can simply define the fields as though they were keyword arguments.
In general, few of the fields are given defaults with the exception of things like trailing edge gap or y-positions, which are always defaulted to zero.
In addition, some parameterization methods have specific values inherent to their methods. These are exposed to the user for convenience, but are also defaulted to the values inherent in the parameterization method.

### Coordinates
All coordinates are given clockwise from the trailing edge, even if the coordinates are split between upper and lower sides.  In other words, the coordinates are given from the lower trailing edge to the leading edge, then the leading edge back to the upper side trailing edge.  In general, if the coordinates are given in an upper and lower split, the leading edge point is repeated.

------------------------------------------------------------------------------------------
# NACA Parameterizations

## NACA 4-series

We begin with the standard NACA 4-series airfoil defined by maximum camber, position of maximum camber, and maximum thickness.

```@example naca4
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

parameters = NACA4(; maximum_camber=2.0, maximum_camber_position=4.0, maximum_thickness=12.0, blunt_trailing_edge=false)
```

!!! note "Defaults"
    The `NACA4` parameter type is the only one with fully defined defaults, which happen to default to the NACA 2412 airfoil with a sharp trailing edge due to its ubiquity.

We can then determine the x,y coordinates from the parameters.

```@example naca4

x, y = naca4(parameters)

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

If you are already familiar with a specific parameterization method, you can also forego defining the parameter object and call the method directly.

```@example naca4
x, y = naca4(2.0, 4.0, 12.0)

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="Direct", linewidth=4, linestyle=:dash)
```

If you have a set of coordinates, and would like to fit NACA 4-series parameters to them, you can use the `determine_naca4` function.

```@example naca4
fit_parameters = determine_naca4(x,y)

x, y = naca4(fit_parameters)

plot!(x, y; label="Fit")
```

```@autodocs; Private=false
Modules = [AirfoilTools]
Pages = ["src/AirfoilTools/parameterizations/naca4.jl"]
```


## NACA 65-series

The NACA 65-series of airfoils was historically used for various linear cascade experiments done by NACA (see ["NACA Report No 824 Summary of Airfoil Data" by  Ira H. Abbott, Albert E. Von Doenhoff, and Louis S. Stivers, Jr. ] (https://ntrs.nasa.gov/citations/19930090976)).
We have implemented some NACA 65-series parameterizations associated with various NACA experimental studies.
In general, the 65-series requires a target lift coefficient, a mean line designation, and a series number.


```@example naca65
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

#define design lift coefficient
clo = 1.0

#define mean line designation
a = 1.0

#define series number
series_number = "3-018"

x, y = AirfoilTools.naca65(clo, a, series_number)

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

There exists a special scaled case for the NACA 65-010 family of airfoils. This special case only requires the user to input a design lift coefficient. using

```@example naca65

#define design lift coefficient
clo = 1.0

x, y = AirfoilTools.naca65_scaled(clo)

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

```@autodocs; Private=false
Modules = [AirfoilTools]
Pages = ["src/AirfoilTools/parameterizations/naca65.jl"]
```

------------------------------------------------------------------------------------------
# Conformal Mapping Parameterizations

Conformally mapped airfoils are defined from analytical solutions for flow about a cylinder.
Currently we have the Joukowsky conformal mapping implemented, which features a cusped trailing edge.  In development is the Karman-Trefftz conformal mapping, which allows for a non-zero trailing edge angle.

The Joukowsky map is based on the center and radius of a circle

```@example joukowsky
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

center = [-0.1; 0.1]
radius = 1.0
pc = plot(
    cosd.(0:360) * radius .+ center[1],
    sind.(0:360) * radius .+ center[2];
    aspectratio=1,
    label="",
    xlabel=L"\xi",
    ylabel=L"\eta",
    framestyle=:origin
)
scatter!(pc, [center[1]], [center[2]]; color=1, label="")

parameters = Joukowsky(; center=center, radius=radius)

x, y = joukowsky(parameters; N=161)

paf = plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")

plot(pc, paf; size=(900, 300))
```

Since conformally mapped airfoils are defined from analytic solutions, we can also use those solutions for validation as done in the [Additional Examples](@ref).

```@example joukowsky
angle_of_attack = 5.0

surface_velocity, surface_pressure_coefficient, cl = joukowsky_flow(
    center, radius, angle_of_attack
)

plot(
    x[2:(end - 1)],
    surface_pressure_coefficient[2:(end - 1)];
    xlabel=L"x",
    ylabel=L"c_p",
    yflip=true,
    label="",
)
```

```@autodocs; Private=false
Modules = [AirfoilTools]
Pages = ["src/AirfoilTools/parameterizations/joukowsky.jl"]
```



------------------------------------------------------------------------------------------
# Class Shape Transformation (CST) Paramterizations

------------------------------------------------------------------------------------------
# Direct Spline Parameterizations

------------------------------------------------------------------------------------------
# PARSEC Parameterizations

------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------

## [C](#)lass [S](#)hape [T](#)ransformations ([CST](#))

The `CST` parameter type is defined as

```@docs
AirfoilTools.CST
```

A CST airfoil can be generated with

```@docs
AirfoilTools.cst
```

For a given set of x,y coordinates a best-fit CST airfoil can be found using

```@docs
AirfoilTools.determine_cst
```

------------------------------------------------------------------------------------------

## Basic B-Spline

The basic B-Spline parameterization comes from ["Universal Airfoil Parametrization Using B-Splines"](https://arc.aiaa.org/doi/10.2514/6.2018-3949) by Rajnarayan, Ning, and Mehr.
It is a cubic B-Spline parameterization based on leading edge radius, trailing edge camber and wedge angle, and optional trailing edge gap distance.

```@docs
AirfoilTools.BasicBSpline
```

A basic B-Spline airfoil can be generated with

```@docs
AirfoilTools.basic_bspline
```

------------------------------------------------------------------------------------------

## [Par](#)ametric [Sec](#)tion ([PARSEC](#))


### Standard
The nominal `PARSEC` parameter type implemented is defined as

```@docs
AirfoilTools.PARSEC
```

A PARSEC airfoil can be generated with

```@docs
AirfoilTools.parsec
```

For a given set of x,y coordinates a best-fit PARSEC airfoil can be found using

```@docs
AirfoilTools.determine_parsec
```

### Modified

Also implemented in AirfoilTools is a modified PARSEC parameterization that give direct control to the trailing edge surfaces of the upper and lower sides.  The `ModifiedPARSEC` type is defined as

```@docs
AirfoilTools.ModifiedPARSEC
```

A Modified PARSEC airfoil can be generated with

```@docs
AirfoilTools.modified_parsec
```

For a given set of x,y coordinates a best-fit Modified PARSEC airfoil can be found using

```@docs
AirfoilTools.determine_modified_parsec
```

## Contributing

We welcome the addition of other common parameterizations.
Additions should have outputs consistent with current parameterizations.
