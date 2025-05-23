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

## Hess-Smith 2D Method

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

## NeuralFoil Lite

[NeuralFoil](https://github.com/peterdsharpe/NeuralFoil) is a multi-layer perceptron model of Xfoil developed by Peter Sharpe in Python.
For a detailed description, see the linked code and his PhD thesis: ["Accelerating Practical Engineering Design Optimization with Computational Graph Transformations"](https://dspace.mit.edu/handle/1721.1/157809?show=full)).
We have ported over a small, basic portion of NeuralFoil--specifically portions providing functionality similar to NeuralFoil's [`get_aero_from_coordinates`](https://github.com/peterdsharpe/NeuralFoil/blob/537860c33496af2fd878999029c59d3272244107/neuralfoil/main.py#L384) function.
(Note that we implement our own least squares solver to determine the Kulfan parameter inputs to the neural net from the input coordinates, thus we have only ported over the neural net portions of neuralfoil and have not included any functionality from [AeroSandbox](https://github.com/peterdsharpe/AeroSandbox) as NeuralFoil does.)
Accessing this NeuralFoil functionality in FLOWFoil can be done using the `NeuralFoil` method:

```@docs
FLOWFoil.NeuralFoil
```

```@example neuralfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

reynolds = 2e6
mach = 0.0

method = NeuralFoil(reynolds, mach; model_size="xlarge")

outputs = analyze([x y], angles_of_attack; method=method)
```

Note that the NeuralFoil method does not allow multi-body analysis like the other methods do as it is based specifically on Xfoil.  We also return a unique output type for the NeuralFoil method that includes all outputs from the NeuralFoil neural net:

```@docs
FLOWFoil.NeuralOutputs
```

!!! note
    The NeuralFoil method here is NOT a replacement for the excellent, and far richer Python implementation.  Rather, this is simply a generally algorithmically differentiable version of the core of NeuralFoil.
    Thus far, we have only tested using ForwardDiff.jl, but we expect other Julia-based automatic differentiation packages to work with little to no modification of the code.
