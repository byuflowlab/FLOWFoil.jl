# Basic Tutorials

<!-- FLOWFoil includes various panel method implementations that are available based on the `method` keyword argument. -->
<!-- Here we go over the available methods and their various options -->

<!-- ## Xfoil Method -->

<!-- An Xfoil-like method, actually based on [mfoil](https://websites.umich.edu/~kfid/codes.html) can be accessed using the `Mfoil` method type: -->

<!-- ```@docs -->
<!-- FLOWFoil.Mfoil -->
<!-- ``` -->

<!-- Note that we have also set `Xfoil=Mfoil` so you can also use the `Xfoil` method type with identical results. -->

<!-- ```@example mfoil -->
<!-- using FLOWFoil -->

<!-- x, y = AirfoilTools.naca4() -->

<!-- angles_of_attack = range(-5.0, 15.0, step=1) -->

<!-- method = Mfoil() -->

<!-- outputs = AirfoilTools.analyze(x, y, angles_of_attack; method=method) -->
<!-- ``` -->

## Lewis' Method for Axisymmetric Bodies

A axisymmetric method based on that described by [Lewis](https://doi.org/10.1017/CBO9780511529542) can be accessed using the `Lewis` method type:

```@docs
FLOWFoil.Lewis
```

```@example lewis
using FLOWFoil

# - DUCT - #

x, r = AirfoilTools.naca4()

# give the duct some diameter
r .+= 1.0

# indicate that the body is not a body of revolution (lying on the axis)
method = Lewis(body_of_revolution = [false])

# angle of attack needs to be zero since it's axisymmetric
outputs = analyze(x, r, 0.0; method=method)
```

!!! note
    No part of the geometry for an axisymmetric body can reside below z=0, otherwise an error will be thrown.


<!-- ## Martensen's Method for Periodic Bodies (Cascades) -->

<!-- A periodic method for cascade analysis based on that developed by [Martensen](https://archive.org/details/nasa_techdoc_19710021012) can be accessed using the `Martensen` method type: -->

<!-- ```@docs -->
<!-- FLOWFoil.Martensen -->
<!-- ``` -->

<!-- ```@example martensen -->
<!-- using FLOWFoil -->

<!-- # - DUCT - # -->

<!-- x, y = AirfoilTools.naca4(6,4,12) -->

<!-- angles_of_attack = range(-5.0, 15.0, step=1) -->

<!-- method = Martensen(pitch=2.0*pi/15.0) -->

<!-- outputs = AirfoilTools.analyze(x, y, angles_of_attack; method=method) -->
<!-- ``` -->

<!-- ## Hess-Smith 2D Method for Educational Purposes -->

<!-- We also have a version of the Hess-Smith method primarily for educational use that can be accessed with the `HessSmith` method type: -->

<!-- ```@docs -->
<!-- FLOWFoil.HessSmith -->
<!-- ``` -->

<!-- ```@example hesssmith -->
<!-- using FLOWFoil -->

<!-- # - DUCT - # -->

<!-- x, y = AirfoilTools.naca4(6,4,12) -->

<!-- angles_of_attack = range(-5.0, 15.0, step=1) -->

<!-- method = HessSmith() -->

<!-- outputs = AirfoilTools.analyze(x, y, angles_of_attack; method=method) -->
<!-- ``` -->
