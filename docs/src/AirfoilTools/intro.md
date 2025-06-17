# Airfoil Tools

AirfoilTools is a submodule of FLOWFoil containing useful airfoil geometry generation and manipulation routines that are commonly used in the BYU FLOW Lab.

## Airfoil Parameterizations

|Method|Status|Fit|
|---|---|---|
|Basic B-Spline|✅|🚧|
|B-Spline with perturbations|🚧|⭕️|
|CST (Kulfan)|✅|✅|
|CST (Circular-arc)|✅|⭕️|
|Joukowsky|✅|⭕️|
|Karman-Trefftz|🚧|⭕️|
|NACA 4-series|✅|✅|
|NACA 65-series|✅|⭕️|
|PARSEC|✅|✅|
|PARSEC (modified)|✅|✅|

\*Fit indicates there is a method implemented to determine parameters from coordinates.

Key:
- ✅ Ready to use
- 🚧 Under Development
- ⭕️ Needs to be added


## Airfoil Geometry Manipulation
- split into upper/lower sides
- re-panel (spline + interpolate)
- flip, rotate, translate, and normalize


## Airfoil Polar Manipulation
- Lift polar linearization


## Contribution

We welcome additional airfoil parameterization methods as well as convenient geometry manipulation routines that make life easier when working with airfoil analysis and optimization.
