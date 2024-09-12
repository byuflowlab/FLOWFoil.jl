# Advanced Examples

## Multiple Airfoil Inviscid Solution

For a multi-element airfoil system, the procedure is identical to a single body system, except we input a vector of matrices for the coordinates of the various bodies.
For this case, we'll use data that comes from ["An Exact Test Case for the Plane Potential Flow About Two Adjacent Lifting Aerofoils" by B. R. Williams.](https://reports.aerade.cranfield.ac.uk/handle/1826.2/2993)

```@example
using FLOWFoil

# SET UP GEOMETRY
af_geom_path = normpath(joinpath(splitdir(pathof(FLOWFoil))[1], "..", "docs", "src", "assets", "two_inviscid_airfoils.jl"))
include(af_geom_path)

# arbitrarily pick an angle of attack
alpha = 0.0

outputs = analyze([[ximain etamain], [xiflap etaflap], alpha)

nothing #hide
```

We see excellent agreement with the analytical solution.

![](assets/two_inviscid_airfoils.jpeg)


## Advanced Viscous Airfoil Options


---

## Axisymmetric Body of Revolution

For this example, we use data from chapter 4 of ["Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems" by  R. I. Lewis](https://doi.org/10.1017/CBO9780511529542)

```@example bor
using FLOWFoil
using PyPlot

include("../../test/data/bodyofrevolutioncoords.jl")

nothing #hide
```

In order to have the correct flags associated with the axisymmetric solver, we generate the mehs using the `generate_axisym_mesh` function.
In addition, since we are modeling a body of revolution, that is, we have an open geometry at the axis of rotation, we need to include the `bodyofrevolution` keyword argument.


![](assets/bodyofrevolution.jpg)

## Axisymmetric Annular Airfoil (Duct)

If we define an airfoil shape in an axisymmetric scheme, we model an annular airfoil, or in other words, a duct.  To do so, we follow a similar procedure to bodies of revolution with the exception that we set `bodyofrevolution=false`.

```@example aa
using FLOWFoil

include("../../test/data/naca_662-015.jl")

nothing #hide
```

As above, we plot experimental results along with our calculated values.

![](assets/annular_airfoil.jpg)


## Axisymmetric Mutli-element Systems

As an example of an multi-element axisymmetric system (such as that used for a ducted rotor), we will simply combine the two previous cases.
We proceed in the same manner for 2D (planar) multi-element systems in that we simply put the various mesh objects in an array together when defining the problem object.

```@example dr
using FLOWFoil

# create annular airfoil mesh object
include("../../test/data/naca_662-015.jl")

# create body of revolution mesh object
include("../../test/data/bodyofrevolutioncoords.jl")

nothing #hide
```

Plotting the geometry and the output velocities and pressures show expected behavior when combining these two cases.

![](assets/ducthubgeom.jpg)

![](assets/multi_body_vs.jpg)

![](assets/multi_body_cp.jpg)
