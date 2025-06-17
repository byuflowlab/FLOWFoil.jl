# Airfoil Tools

AirfoilTools is a submodule of FLOWFoil containing useful airfoil geometry generation and manipulation routines that are commonly used in the BYU FLOW Lab.

## Airfoil Generation

Several common airfoil parameterization methods are implemented in AirfoilTools for generating various airfoil geometries.

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

Occasionally, we might want to manipulate airfoil geoemtries for some reason.  We also have methods for the following in the AirfoilTools module:

- split into upper/lower sides
- flip, rotate, translate, and normalize
- re-panel (spline + interpolate)


## Contribution

We welcome additional airfoil parameterization methods as well as convenient geometry manipulation routines that make life easier when working with airfoil analysis and optimization.
