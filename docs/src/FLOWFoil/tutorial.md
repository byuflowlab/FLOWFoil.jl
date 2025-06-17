# Quick Start

Running FLOWFoil can be done simply with a single method: `analyze`
As an introductory example, we will do a quick analysis of a NACA 2412 airfoil, with coordinates from one of the available methods in the [Airfoil Tools](@ref) sub-module.

```@example tutorial
using FLOWFoil

# note: 2412 is default for the NACA 4-series implemenation
x, y = AirfoilTools.naca4()
nothing #hide
```

Let's plot the geometry.
```@example tutorial
# plot geometry
include("../assets/plots_default.jl") #hide
using Plots
using LaTeXStrings
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

!!! note
    AirfoilTools generates airfoil coordinates in the format required for FLOWFoil: starting at the trailing edge, and proceeding clockwise around the airfoil.


Let's finalize the required inputs and run the analysis.
```@example tutorial
# choose one or more angles of attack
angles_of_attack = range(-5.0, 15.0, step=1)

# analyze
outputs = analyze(x, y, angles_of_attack)
```

And then we can plot some outputs, for example, the pressure distribution.
```@example tutorial
# plot pressure distribution at the 6th angle of attack.
plot(x, outputs.cp[:, 6]; xlabel=L"x", ylabel=L"c_p", yflip=true, label="")
```

Or the lift and drag polars.
```@example tutorial
# plot lift
pcl = plot(angles_of_attack, outputs.cl, xlabel=L"\alpha", ylabel=L"c_\ell", label="")
pcd = plot(angles_of_attack, outputs.cd, xlabel=L"\alpha", ylabel=L"c_d", label="")
plot(pcl, pcd; size=(900,300))
```

Output structures depend on the method selected, but in general you'll get lift (cl) and drag (cd) and some other values depending on the method.
In this case, we are using the default Xfoil/[Mfoil](@ref) method and the outputs are an `InviscidOutputs` type, which all of the non-wrapped methods use currently:

```@docs
FLOWFoil.InviscidOutputs
```

See [Tutorials](@ref) for the various methods and additional output structure types.
