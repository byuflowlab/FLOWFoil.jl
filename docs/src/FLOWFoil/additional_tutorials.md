# Tutorials

FLOWFoil includes various panel method implementations that are available based on the `method` keyword argument.
The analyze function is the way to run analyses with any method available in FLOWFoil, you just need to select the method you want and provide any additionally required inputs.

```@docs
FLOWFoil.analyze
```

## Xfoil Method

An Xfoil-like method, based on [Xfoil](https://websites.umich.edu/~kfid/codes.html), can be accessed using the `Xfoil` method type. It is also exported as `Mfoil` (an alias) since the implementation follows Drela's `mfoil` derivation. The method has two modes:

- **Inviscid** (default): a linear-vortex panel solve over the airfoil. Returns an `InviscidOutputs` object.
- **Viscous**: a coupled boundary-layer + inviscid solve. Set `viscous=true` and provide `Re`, `Mach`, and the other viscous knobs. Returns a NamedTuple (see [Viscous Mode](@ref) below).

### Inviscid Mode

```@example Xfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

method = Xfoil()

outputs = analyze(x, y, angles_of_attack; method=method)
```

The inviscid Xfoil outputs are of type `InviscidOutputs`. This is also the default method used in the [Quick Start](@ref).

### Viscous Mode

When `viscous=true` is set, `analyze` runs the coupled mfoil-style boundary-layer + inviscid panel solve described by [Drela](https://doi.org/10.1007/978-3-642-84010-4_1) and [Fidkowski](https://doi.org/10.2514/1.J061341). Each angle of attack is solved independently via an internal loop. A worked example with plots lives on the [Additional Examples](@ref) page.

```@example xfoil_viscous
using FLOWFoil

x, y = AirfoilTools.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)

method = Xfoil(viscous=true, Re=1e6, Mach=0.0, ncrit=9.0, etol=1e-10, maxiters=50)

outputs = analyze([x y], [5.0]; method=method)
nothing # hide
```

```@example xfoil_viscous
(cl=outputs.cl[1],
 cd=outputs.cd[1],
 cm=outputs.cm[1],
 transition_upper=outputs.transition_upper[1],
 transition_lower=outputs.transition_lower[1],
 converged=outputs.converged[1])
```

The viscous return is a flat NamedTuple. The top-level scalar fields per angle of attack are `cl`, `cd`, `cm`, `transition_upper`, `transition_lower`, and `converged`. The inviscid intermediates `gamma_ref` (node vortex strengths) and `A` (the linear-system influence matrix) are also exposed alongside. Per-station boundary-layer state for each α is in `outputs.per_alpha[i]` and includes `thetas_u/l/w`, `deltas_u/l/w`, `state3s_u/l/w` (laminar amplification factor `n` or turbulent shear-stress coefficient `cτ`), and the edge velocities `Veu`, `Vel`, `Vwake`.

!!! note "Viscous-mode constraints"
    - **Single body only.** A viscous run on a multi-body system raises an error.
    - **Blunt trailing edge required.** Use `AirfoilTools.naca4(...; blunt_trailing_edge=true)` or supply your own blunt-TE coordinates.

```@docs
FLOWFoil.Xfoil
```

## Lewis' Method for Axisymmetric Bodies

An axisymmetric method based on the one described by [Lewis](https://doi.org/10.1017/CBO9780511529542) can be accessed using the `Lewis` method type:

```@example lewis
using FLOWFoil

# Since this is an axisymmetric method, we'll use r instead of y
x, r = AirfoilTools.naca4()

# give the duct some diameter so it doesn't have negative radial dimensions (see warning below)
r .+= 1.0

# indicate that the body is not a body of revolution (i.e. a duct)
method = Lewis(; body_of_revolution=false)

# note: we need to input an an angle of attack, even though it is unused in this method.
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
    If the `cascade` option is set to false, this method becomes a standard planar airfoil method, but uses constant vortices, so the Xfoil method is the superior method in that case.

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
