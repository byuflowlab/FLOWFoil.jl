# FLOWFoil.jl ([Fl]()ight, [O]()ptimization, and [W]()ind Air[Foil]() Analysis)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://byuflowlab.github.io/FLOWFoil.jl/stable)
[![Build Status](https://github.com/byuflowlab/FLOWFoil.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/FLOWFoil.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


FLOWFoil is a collection of two dimensional potential flow solvers (panel methods) for airfoils, airfoil systems, and axisymmetric sections/systems.

The following table includes the list of currently available methods (usage can be found in the docs).

|Method|Inviscid Mfoil/Xfoil|Lewis|Martensen|LegacyXfoil|NeuralFoil|
|---|---|---|---|---|---|
Description|A re-derivation of the Mfoil/Xfoil method and implementation in Julia. Currently only the inviscid portions are derived/implemented. | An axisymmetric method, useful for ducts and bodies of revolution. | A periodic method (with optional planar functionality) for axial cascades. |  Wrapper of Xfoil.jl | Wrapper of NeuralFoil.jl |
Panel Type | Linear Vortex (+ Constant Source) | Constant Vortex | Constant Vortex | Constant Source + Single Vortex | Linear Vortex + Constant Source | Multi-Layer Perceptron |
Boundary Condition|Dirichlet|Dirichlet|Dirichlet|Dirichlet|N/A|
Viscous|🚧|⭕️|⭕️|✅|✅|
Single-body Functional|✅|✅|✅|✅|✅|
Multi-body Functional|🚧|✅|🚧|❌|❌|
Able to model blunt trailing edges|✅|⭕️|⭕️|✅|✅|
AD Compatible (ForwardDiff.jl)|✅|✅|✅|❌|✅|
References or Wrapped Packages|[1](https://websites.umich.edu/~kfid/codes.html), [2](https://web.mit.edu/drela/Public/papers/xfoil_sv.pdf)|[3](https://doi.org/10.1017/CBO9780511529542) |[3](https://doi.org/10.1017/CBO9780511529542) | [Xfoil.jl](https://github.com/byuflowlab/Xfoil.jl)| [NeuralFoil.jl](https://github.com/byuflowlab/NeuralFoil.jl)  |

Key:
- ✅ Implemented
- 🚧 Under Development
- ⭕️ Will likely not be implemented
- ❌ Cannot be implemented

## References:

1. [Fidkowski, K. J., "A Coupled Inviscid-Viscous Airfoil Analysis Solver, Revisited," AIAA Journal, 2021.](https://doi.org/10.2514/1.J061341)
2. [Drela, M., "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils," 1989.](https://doi.org/10.1007/978-3-642-84010-4_1)
3. [R. I. Lewis, "Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems," 1991](https://doi.org/10.1017/CBO9780511529542)
