# Quick Start

Running FLOWFoil can be done simply with a single method: `analyze`

```@docs
FLOWFoil.analyze
```

As an introductory example, we will do a quick analysis of a NACA 2412 airfoil, with coordinates from one of the available methods in the [Airfoil Tools](@ref) sub-module.

!!! note
    For any airfoil coordinate generation method (from FLOWFoil or otherwise), the coordinates must start at the trailing edge, and proceed clockwise around the airfoil.

```@julia
using FLOWFoil

# 2412 is default
x, y = AirfoilTools.naca4()

# choose one or more angles of attack
angles_of_attack = range(-5.0, 15.0, step=1)

outputs = analyze(x, y, angles_of_attack)
```

The `outputs` named tuple will contain
- `cl`: lift coefficient of each body
- `cd`: total drag coefficient of each body
- `cdp`: profile drag coefficient of each body
- `cm`: moment coefficient of each body
- `tangential_velocities`: surface velocities on each body
- `surface_pressures`: surface pressures on each body
- `convergenced`: convergence flag
- `auxiliary outputs`: a named tuple that contains additional outputs applicable to the method used.
