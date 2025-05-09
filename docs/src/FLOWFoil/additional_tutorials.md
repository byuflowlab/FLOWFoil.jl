# Basic Tutorials

FLOWFoil includes various panel method implementations that are available based on the `method` keyword argument.
Here we go over the available methods and their various options

## Xfoil Method

An Xfoil-like method, actually based on [mfoil](https://websites.umich.edu/~kfid/codes.html) can be accessed using the `Mfoil` method type:

```@docs
FLOWFoil.Mfoil
```

Note that we have also set `Xfoil=Mfoil` so you can also use the `Xfoil` method type with identical results.
Currently, this method only includes the inviscid parts of Xfoil/Mfoil, so Reynolds and Mach number inputs do nothing if used.

```@example Mfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

# viscous solver not yet implemented
method = Mfoil(viscous=false)

outputs = analyze(x, y, angles_of_attack; method=method)
```

## Lewis' Method for Axisymmetric Bodies

An axisymmetric method based on that described by [Lewis](https://doi.org/10.1017/CBO9780511529542) can be accessed using the `Lewis` method type:

```@docs
FLOWFoil.Lewis
```

```@example lewis
using FLOWFoil

x, r = AirfoilTools.naca4()

# give the duct some diameter (see note below)
r .+= 1.0

# indicate that the body is not a body of revolution (i.e. a duct)
method = Lewis(; body_of_revolution=[false])

# note: we need to input an an angle of attack, even though it is unused
outputs = analyze(x, r, [0.0]; method=method)
```

!!! note
    No part of the geometry for an axisymmetric body can reside below z=0, otherwise an error will be thrown.

## Martensen's Method for Periodic Bodies (Cascades)

A periodic method for cascade analysis based on that developed by [Martensen](https://archive.org/details/nasa_techdoc_19710021012) can be accessed using the `Martensen` method type:

```@docs
FLOWFoil.Martensen
```
!!! note
    If the `cascade` option is set to false, this method becomes a standard planar airfoil method, but uses constant rather than linear vortices, so the Mfoil/Xfoil method is superior in that case.

```@example martensen
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

# The cascade method requires solidity (closeness) of sections and stagger (inflow angle - angle of attack)
method = Martensen(solidity=1.2, stagger=15.0)

outputs = analyze(x, y, angles_of_attack; method=method)
```

## Hess-Smith 2D Method for Educational Purposes

We also have a version of the Hess-Smith method primarily for educational use that can be accessed with the `HessSmith` method type:

```@docs
FLOWFoil.HessSmith
```

```@example HessSmith
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

method = HessSmith(V_inf=1.0)

outputs = analyze(x, y, angles_of_attack; method=method)
```
