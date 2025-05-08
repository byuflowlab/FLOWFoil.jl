# Basic Tutorials

FLOWFoil includes various panel method implementations that are available based on the `method` keyword argument.
Here we go over the available methods and their various options

## Xfoil Method

An Xfoil-like method, actually based on [mfoil](https://websites.umich.edu/~kfid/codes.html) can be accessed using the `Mfoil` method type:

```@docs
FLOWFoil.Mfoil
```

Note that we have also set `Xfoil=Mfoil` so you can also use the `Xfoil` method type with identical results.
Currently, this method only includes the inviscid parts of Xfoil/Mfoil.

```@example Mfoil
using FLOWFoil

#setup x and y coordinates for the airfoil
x, y = AirfoilTools.naca4()

#create a range of angles of attack
angles_of_attack = range(-5.0, 15.0, step=1)

#solve for invsicid case
viscous = false

#set up Mfoil method
method = Mfoil(viscous)

#solve for outputs
outputs = analyze(x, y, angles_of_attack; method=method)
```

## Lewis' Method for Axisymmetric Bodies

A axisymmetric method based on that described by [Lewis](https://doi.org/10.1017/CBO9780511529542) can be accessed using the `Lewis` method type:

```@docs
FLOWFoil.Lewis
```

```@example lewis
using FLOWFoil

#setup airfoil coordinates
x, r = AirfoilTools.naca4()

# give the duct some diameter
r .+= 1.0

# indicate that the body is not a body of revolution (lying on the axis)
method = Lewis(body_of_revolution = [false])

# angle of attack defaults to zero, which is what we want for the axisymmetric case
outputs = analyze(x, r; method=method)
```

!!! note
    No part of the geometry for an axisymmetric body can reside below z=0, otherwise an error will be thrown.


## Martensen's Method for Periodic Bodies (Cascades)

A periodic method for cascade analysis based on that developed by [Martensen](https://archive.org/details/nasa_techdoc_19710021012) can be accessed using the `Martensen` method type:

```@docs
FLOWFoil.Martensen
```

```@example martensen
using FLOWFoil

#setup airfoil coordinates
x, y = AirfoilTools.naca4(6,4,12)

#specify range of angle of attack
inflow_angles = range(-5.0, 15.0, step=1)

#solve for airfoil cascade
cascade = true

#setup solidity - only important in cascade case. If cascade = false, then value of solidity will not effect the outputs
solidity = 1.0

#setup airfoil stagger (inflow angle minus angle of attack). When in doubt set this to 0.0
stagger = 0.0

#Value of solidity that the cascade method will transition to solving using the planar method (ie single airfoil instead of a cascade). In this case it is used for blending the cascade solution and planar solution together when the solidity is close to the transition value. When in doubt set this to 1e-4.
transition_value = 0.0001

#curvature correction is currently not working, always set to false
curvature_correction = false

#specify Martensen method
method = Martensen(cascade,solidity,stagger, transition_value, curvature_correction)

#solve for outputs
outputs = analyze(x, y, inflow_angles; method=method)
```

## Hess-Smith 2D Method for Educational Purposes

We also have a version of the Hess-Smith method primarily for educational use that can be accessed with the `HessSmith` method type:

```@docs
FLOWFoil.HessSmith
```

```@example HessSmith
using FLOWFoil

#specify airfoil coordinates
x, y = AirfoilTools.naca4(6,4,12)

#specify range of angle of attack
angles_of_attack = range(-5.0, 15.0, step=1)

#specify magnitude of the free stream velocity
vinf = 1.0

#specify HessSmith method
method = HessSmith(vinf)

#solve for outputs
outputs = analyze(x, y, angles_of_attack; method=method)
```
