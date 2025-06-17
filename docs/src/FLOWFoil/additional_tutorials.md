# Tutorials

FLOWFoil includes various panel method implementations that are available based on the `method` keyword argument.
The analyze function is the way to run analyses with any method available in FLOWFoil, you just need to select the method you want and provide any additionally required inputs.

```@docs
FLOWFoil.analyze
```

## Mfoil (inviscid) Method

An Xfoil-like method, actually based on [mfoil](https://websites.umich.edu/~kfid/codes.html) can be accessed using the `Mfoil` method type:

Note that we have also set `Xfoil=Mfoil` so you can also use the `Xfoil` method type with identical results.
Currently, this method only includes the inviscid methods of Xfoil/Mfoil.

```@example Mfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

# viscous solver not yet implemented
method = Mfoil(viscous=false)

outputs = analyze(x, y, angles_of_attack; method=method)
```

Currently, the Mfoil method outputs are of type InviscidOutputs.  This is also the default method used in the [Quick Start](@ref).

```@docs
FLOWFoil.Mfoil
```

## Lewis' Method for Axisymmetric Bodies

An axisymmetric method based on that described by [Lewis](https://doi.org/10.1017/CBO9780511529542) can be accessed using the `Lewis` method type:

```@example lewis
using FLOWFoil

# Since this is an axisymmetric method, we'll use r instead of y
x, r = AirfoilTools.naca4()

# give the duct some diameter so it doesn't have negative radial dimensions (see warning below)
r .+= 1.0

# indicate that the body is not a body of revolution (i.e. a duct)
method = Lewis(; body_of_revolution=false)

# note: we need to input an an angle of attack, even though it is unused
outputs = analyze(x, r, 0.0; method=method)
```

The comments here mention multiple bodies, for more information, see the multi-body example: [Axisymmetric Mutli-element Systems](@ref) on the next page.

The outputs for the Lewis method are also of type InviscidOutputs.

```@docs
FLOWFoil.Lewis
```

!!! warning
    No part of the geometry for an axisymmetric body can reside below z=0, otherwise an error will be thrown.

## Martensen's Method for Axial Cascades

A periodic method for cascade analysis based on that developed by [Martensen](https://archive.org/details/nasa_techdoc_19710021012) can be accessed using the `Martensen` method type:

```@example martensen
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

# The cascade method requires solidity (closeness) of sections and stagger (inflow angle - angle of attack)
method = Martensen(solidity=1.2, stagger=15.0)

outputs = analyze(x, y, angles_of_attack; method=method)
```

The InviscidOutputs type is also used for the Martensen method.

```@docs
FLOWFoil.Martensen
```
!!! note
    If the `cascade` option is set to false, this method becomes a standard planar airfoil method, but uses constant vortices, so the Mfoil/Xfoil method is the superior method in that case.

## NeuralFoil Method

[NeuralFoil](https://github.com/peterdsharpe/NeuralFoil) is a multi-layer perceptron model of Xfoil.
We provide the Neuralfoil Method through the `NeuralFoil` method type and is accessed through the [NeuralFoil.jl](https://github.com/byuflowlab/NeuralFoil.jl) package:

```@example neuralfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

reynolds = 2e6
mach = 0.0

method = NeuralFoil(reynolds, mach; model_size="xlarge", n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)

outputs = analyze([x y], angles_of_attack; method=method)
```

```@docs
FLOWFoil.NeuralFoil
FLOWFoil.NeuralFoil(reynolds)
```

Note that the NeuralFoil method does not allow multi-body analysis like the other methods do as it is based specifically on Xfoil.  We also return a separate output type for the NeuralFoil method from the NeuralFoil.jl namespace.

## LegacyXfoil Method

We also have the LegacyXfoil method that is based on Xfoil and can be accessed with the `LegacyXfoil` method type:

```@example legacyxfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

reynolds = 2e6

method = LegacyXfoil(reynolds; npan=140)

outputs = analyze([x y], angles_of_attack; method=method)
```

```@docs
FLOWFoil.LegacyXfoil
FLOWFoil.LegacyXfoil(reynolds)
```

Note that we return a separate output type for the LegacyXFoil method:

```@docs
FLOWFoil.LegacyXFOutputs
```
