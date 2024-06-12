# Common Parameterizations

Included in the AirfoilTools module are the following common airfoil parameterizations for convenience in generating various airfoil geometries.

### Parameter Types
Each of the parameterizations below have associated composite types whose fields are the parameters used in defining the airfoil geometries.  Each of these composite types is defined using the `@kwdef` macro such that the user does not need to remember the order of the fields, but can simply define the fields as though they were keyword arguments.
In general, few of the fields are given defaults with the exception of thigns like trailing edge gap or z-positions, which are always defaulted to zero.
In addition, some parameterization methods have specific values inherent to their methods. These are exposed to the user for convenience, but are also defaulted to the values inherent in the parameterization method.

### Coordinates
For all of the parameterizations, the x-coordinates are defined to be along the chord-wise direction, and the z-coordinates are orthogonal to the x-coordinates.  This is done to be similar to the standard airframe coordinate system, with the x-direction being positve toward the back of the airframe, and the z-direction being up.

In addition all coordinates are given clockwise from the trailing edge, even if the coordinates are split between upper and lower sides.  In other words, the coordinates are given from the lower trailing edge to the leading edge, then the leading edge back to the upper side trailing edge.  In general, if the coordinates are given in an upper and lower split, the leading edge point is repeated.

------------------------------------------------------------------------------------------

## NACA 4-series

We begin with the standard NACA 4-series airfoil, defining the parameters in the `NACA4` type as

```@docs
AirfoilTools.NACA4
```

!!! note "Defaults"
    The `NACA4` parameter type is the only one with fully defined defaults, which happen to default to the NACA 2412 airfoil with a sharp trailing edge due to its ubiquity.

A NACA 4-series airfoil can be defined with either of the `naca4` methods.

```@docs
AirfoilTools.naca4
```

------------------------------------------------------------------------------------------

## Conformal Mapping

AirfoilTools includes two conformal mapping methods; the first is the Joukowsky method, and the second is the KarmanTrefftz method (which is a variation on Joukowsky allowing for non-cusped trailing edges)

### Joukowsky

The `Joukowsky` parameter type is defined as

```@docs
AirfoilTools.Joukowsky
```

A Joukowsky airfoil can be generated with

```@docs
AirfoilTools.joukowsky
```

### Karman-Trefftz

The `KarmanTrefftz` parameter type is defined as

```@docs
AirfoilTools.KarmanTrefftz
AirfoilTools.KarmanTrefftz(center,wedge_angle)
```

A Karman-Trefftz airfoil can be generated with

```@docs
AirfoilTools.karman_trefftz
```

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

### Modified

Also implemented in AirfoilTools is a modified PARSEC parameterization that give direct control to the trailing edge surfaces of the upper and lower sides.  The `ModifiedPARSEC` type is defined as

```@docs
AirfoilTools.ModifiedPARSEC
```

A Modified PARSEC airfoil can be generated with

```@docs
AirfoilTools.modified_parsec
```

## Contributing

We welcome the addition of other common parameterizations.
Additions should have outputs consistent with current parameterizations.
