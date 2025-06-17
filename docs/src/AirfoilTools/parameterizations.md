# Common Parameterizations

### Parameter Types
Most of the parameterizations below have associated composite types whose fields are the parameters used in defining the airfoil geometries.  Each of these composite types is defined using the `@kwdef` macro such that the user does not need to remember the order of the fields, but can simply define the fields as though they were keyword arguments.
In general, few of the fields are given defaults with the exception of things like trailing edge gap or y-positions, which are always defaulted to zero.
In addition, some parameterization methods have specific values inherent to their methods. These are exposed to the user for convenience, but are also defaulted to the values inherent in the parameterization method.

### Coordinates
All coordinates are given clockwise from the trailing edge, even if the coordinates are split between upper and lower sides.  In other words, the coordinates are given from the lower trailing edge to the leading edge, then the leading edge back to the upper side trailing edge.  In general, if the coordinates are given in an upper and lower split, the leading edge point is repeated.
Note that all coordinates are given with x being in the chord-wise direction and y being orthogonal to x (just like a standard Cartesian x,y plot). Airfoil geometries are all given normalized to the chord length, so the x values will generally have the leading edge at x=0 and the trailing edge at x=1, though there may be some slight variation depending on camber and trailing edge gap.

------------------------------------------------------------------------------------------
# NACA Parameterizations

## NACA 4-series

We begin with the standard [NACA 4-series airfoil](https://en.wikipedia.org/wiki/NACA_airfoil) defined by maximum camber, position of maximum camber, and maximum thickness.

```@example naca4
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

parameters = NACA4(;
    maximum_camber=2.0,
    maximum_camber_position=4.0,
    maximum_thickness=12.0,
    blunt_trailing_edge=false,
)
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

plot(
    x,
    y;
    aspectratio=1,
    xlabel=L"x",
    ylabel=L"y",
    label="Direct",
    linewidth=4,
    linestyle=:dash,
)
```

If you have a set of coordinates, and would like to fit NACA 4-series parameters to them, you can use the `determine_naca4` function.

```@example naca4
fit_parameters = determine_naca4(x,y)

x, y = naca4(fit_parameters)

plot!(x, y; label="Fit")
```

```@docs
AirfoilTools.NACA4
AirfoilTools.naca4
AirfoilTools.determine_naca4
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

```@docs
AirfoilTools.naca65
AirfoilTools.naca65_scaled
```

------------------------------------------------------------------------------------------
# Conformal Mapping Parameterizations

Conformally mapped airfoils are defined from analytical solutions for flow about a cylinder.
Currently we have the [Joukowsky conformal mapping](https://en.wikipedia.org/wiki/Joukowsky_transform) implemented, which features a cusped trailing edge.  In development is the Karman-Trefftz conformal mapping, which allows for a non-zero trailing edge angle.

## Joukowsky

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

parameters = Joukowsky(; center, radius)

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

```@docs
AirfoilTools.Joukowsky
AirfoilTools.joukowsky
AirfoilTools.joukowsky_flow
```

------------------------------------------------------------------------------------------
# Class Shape Transformation (CST) Paramterizations

## Standard Kulfan CST

We implement the CST parameterization presented by [Kulfan](https://doi.org/10.2514/6.2007-62), where airfoil shapes are defined based on upper and lower coefficients.

```@example cst
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

parameters = CST(;
    upper_coefficients=[0.2; 0.3; 0.2; 0.2], lower_coefficients=[-0.1; 0.1; 0.0; 0.0]
)

x, y = cst(parameters)

plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

This CST method also has a fitting function.

```@example cst
x, y = naca4()
plot(
    x,
    y;
    aspectratio=1,
    xlabel=L"x",
    ylabel=L"y",
    label="Coordinates",
    linewidth=4,
    linestyle=:dash,
)

fit_parameters = determine_cst(x, y; n_upper_coefficients=4, n_lower_coefficients=4)

x, y = cst(fit_parameters)
plot!(x, y; label="Fit")
```

```@docs
AirfoilTools.CST
AirfoilTools.cst
AirfoilTools.determine_cst
```

## Circular Arc Camber CST

We also implement a CST-based circular arc camber line airfoil parameterization that may be useful for axial cascade geometry.  This parameterization comes from ["Aerodynamics of Low Reynolds Number Axial Compressor Sections"](https://doi.org/10.2514/6.2015-1934) by Maffioli et al., and is defined from a camber, maximum thickness, and position of maximum thickness.

```@example circle_cst
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

parameters = CircularArcCST(;
    maximum_camber=0.05, maximum_thickness_postition = 0.3, maximum_thickness = 0.12
)

x, y = circular_arc_cst(parameters)
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

------------------------------------------------------------------------------------------
# Spline Parameterizations

## Basic B-Spline

The basic B-Spline parameterization comes from ["Universal Airfoil Parametrization Using B-Splines"](https://arc.aiaa.org/doi/10.2514/6.2018-3949) by Rajnarayan, Ning, and Mehr.
It is a cubic B-Spline parameterization based on leading edge radius, trailing edge camber and wedge angle, and optional trailing edge gap distance.

```@example bspline
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

parameters = BasicBSpline(;
    leading_edge_radius=0.015, trailing_edge_camber_angle=12.0, wedge_angle=12.0
)

x, y = basic_bspline(parameters)
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

```@docs
AirfoilTools.BasicBSpline
AirfoilTools.basic_bspline
```
------------------------------------------------------------------------------------------
# PARametric SECtion (PARSEC) Parameterizations

## Standard PARSEC

We implement the PARSEC method developed by Sobieczky (using [this](https://doi.org/10.1016/j.ast.2013.11.006) reference).
There are 11 parameters in the traditional PARSEC parameterization method. They are as follows:

| Parameter | Definition |
|------|:------|
| $r_{LE}$ | Leading Edge Radius |
| $x_u$ | Chordwise position of maximum thickness of upper surface |
| $x_l$ | Chordwise position of maximum thickness of lower surface |
| $y_u$ | y-coordinate at maximum thickness of upper surface |
| $y_l$ | y-coordinate at maximum thickness of lower surface |
| $y_{xx_u}$ | Second derivative of upper surface curve at point of maximum thickness |
| $y_{xx_l}$ | Second derivative of lower surface curve at point of maximum thickness |
| $\alpha_{TE}$ | Trailing edge angle |
| $\beta_{TE}$ | Boat-tail angle |
| $y_{TE}$ | y-position of center of trailing edge |
| $\Delta y_{TE}$ | y-distance between upper and lower surface trailing edge points |

```@example parsec
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

parameters = PARSEC(;
    leading_edge_radius=0.015,
    maximum_thickness_xu=0.33,
    maximum_thickness_xl=0.20,
    maximum_thickness_yu=0.08,
    maximum_thickness_yl=-0.04,
    curvature_u=-0.63,
    curvature_l=0.30,
    trailing_edge_angle=-0.05,
    boattail_angle=-0.15,
    trailing_edge_y=0.0,
    trailing_edge_gap=0.0,
)

x, y = parsec(parameters)
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

The PARSEC method also has a fit functionality implemented.

```@example parsec
x, y = naca4()
plot(
    x,
    y;
    aspectratio=1,
    xlabel=L"x",
    ylabel=L"y",
    label="Coordinates",
    linewidth=4,
    linestyle=:dash,
)

fit_parameters = determine_parsec(x, y)

x, y = parsec(fit_parameters)
plot!(x, y; label="Fit")
```

```@docs
AirfoilTools.PARSEC
AirfoilTools.parsec
AirfoilTools.determine_parsec
```

## Modified PARSEC

Also implemented in AirfoilTools is a modified PARSEC parameterization that gives direct control to the trailing edge surfaces of the upper and lower sides.
In order to create a more intuitive parameterization, the trailing edge parameters, both position and angles are replaced to be directly the trailing edge positions of the upper and lower surface and trailing edge surface angles.

| Original Parameter | Modified Parameter | New Definition |
|------|------|:------|
| $\alpha_{TE}$ | $\theta_{TE_u}$ | Upper surface trailing edge tangent angle |
| $\beta_{TE}$ | $\theta_{TE_l}$ | Upper surface trailing edge tangent angle |
| $y_{TE}$ | $y_{TE_u}$ | y-coordinate of upper surface trailing edge |
| $\Delta y_{TE}$ | $y_{TE_l}$ | y-coordinate of lower surface trailing edge |

```@example modparsec
using FLOWFoil.AirfoilTools
using Plots
using LaTeXStrings
include("../assets/plots_default.jl") #hide

parameters = ModifiedPARSEC(;
    leading_edge_radius=0.015,
    maximum_thickness_xu=0.33,
    maximum_thickness_xl=0.20,
    maximum_thickness_yu=0.08,
    maximum_thickness_yl=-0.04,
    curvature_u=-0.63,
    curvature_l=0.30,
    trailing_edge_tangent_u=-0.2,
    trailing_edge_tangent_l=0.1,
    trailing_edge_yu=0.0,
    trailing_edge_yl=0.0,
)

x, y = modified_parsec(parameters)
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

The Modified PARSEC method also has a fit functionality implemented.

```@example parsec
x, y = naca4()
plot(
    x,
    y;
    aspectratio=1,
    xlabel=L"x",
    ylabel=L"y",
    label="Coordinates",
    linewidth=4,
    linestyle=:dash,
)

fit_parameters = determine_modified_parsec(x, y)

x, y = modified_parsec(fit_parameters)
plot!(x, y; label="Fit")
```

```@docs
AirfoilTools.ModifiedPARSEC
AirfoilTools.modified_parsec
AirfoilTools.determine_modified_parsec
```

------------------------------------------------------------------------------------------

# Contributing

We welcome the addition of other common parameterizations, as well as improvements/additions to currently implemented options.
Additions should have outputs consistent with current parameterizations.
