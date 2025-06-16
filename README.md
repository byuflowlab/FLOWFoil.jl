# FLOWFoil.jl ([Fl]()ight, [O]()ptimization, and [W]()ind Air[Foil]() Analysis)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://byuflowlab.github.io/FLOWFoil.jl/stable)
[![Build Status](https://github.com/byuflowlab/FLOWFoil.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/FLOWFoil.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


FLOWFoil is a collection of two dimensional potential flow solvers (panel methods) for airfoils, airfoil systems, and axisymmetric sections/systems.
The following table includes the list of currently available methods (usage can be found in the docs).

| Method | Type | Multi-body Compatible | AD Compatible (ForwardDiff.jl) | References or Wrapped Pacakge |
|---|---|---|---|---|
| Mfoil | Inviscid, linear vortex | ✅ | ✅ | [1](https://websites.umich.edu/~kfid/codes.html), [2](https://web.mit.edu/drela/Public/papers/xfoil_sv.pdf) |
| Lewis | Inviscid, axisymmetric, constant vortex | ✅ | ✅ | [3](https://doi.org/10.1017/CBO9780511529542) |
| Martensen | Inviscid, periodic, constant vortex | ✅ | ✅ | [3](https://doi.org/10.1017/CBO9780511529542) |
| HessSmith | Inviscid, constant source, single vortex | ✅ | ✅ | [4](https://byu.box.com/shared/static/ywfayozbj3sr2ot6b32u8nqk5brqvurt.pdf) |
| LegacyXfoil | Wrapper for Xfoil.jl | ❌ | ❌ | [Xfoil.jl](https://github.com/byuflowlab/Xfoil.jl) |
| NeuralFoil | Wrapper for NeuralFoil.jl | ❌ | ✅ | [NeuralFoil.jl](https://github.com/byuflowlab/NeuralFoil.jl)  |

## References:

1. [Fidkowski, K. J., "A Coupled Inviscid-Viscous Airfoil Analysis Solver, Revisited," AIAA Journal, 2021.](https://doi.org/10.2514/1.J061341)
2. [Drela, M., "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils," 1989.](https://doi.org/10.1007/978-3-642-84010-4_1)
3. [R. I. Lewis, "Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems," 1991](https://doi.org/10.1017/CBO9780511529542)
4. [Ning, A. "Computational Aerodynamics," 2022](https://byu.box.com/shared/static/ywfayozbj3sr2ot6b32u8nqk5brqvurt.pdf)

