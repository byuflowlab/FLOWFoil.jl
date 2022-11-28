# FLOWFoil.jl ([Fl]()ight, [O]()ptimization, and [W]()ind Air[Foil]() Analysis)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://byuflowlab.github.io/FLOWFoil.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://byuflowlab.github.io/FLOWFoil.jl/dev)
[![Build Status](https://github.com/byuflowlab/FLOWFoil.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/FLOWFoil.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


FLOWFoil is a two dimensional potential flow solver (panel method) for airfoils, airfoil systems, and axisymmetric sections/systems.

The formulation used for the planar (airfoils) systems in FLOWFoil follows closely those of [XFoil](https://web.mit.edu/drela/Public/web/xfoil/) and [mfoil](http://www-personal.umich.edu/~kfid/codes.html) (see also [references](#References) below).

For axisymmetric sections (bodies of revolution, ducts) and systems, FLOWFoil follows closely the formulations laid out in [Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems](https://doi.org/10.1017/CBO9780511529542) (which also appear to be similar to the methods applied in the [Ducted Fan Design Code](http://web.mit.edu/drela/Public/web/dfdc/), a program used for analysis of low Reynolds number ducted rotors).

Currently, FLOWFoil only has inviscid capabilties for single and multi-element systems.
Future additions will include visous capabilties for both single and multi-element airfoils as well.



## References:

 - [Drela, M., “XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils,” 1989.](https://doi.org/10.1007/978-3-642-84010-4_1)
 - [Fidkowski, K. J., “A Coupled Inviscid-Viscous Airfoil Analysis Solver, Revisited,” AIAA Journal, 2021.](https://doi.org/10.2514/1.J061341)
 - [R. I. Lewis, "Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems," 1991](https://doi.org/10.1017/CBO9780511529542)


$\checkmark$
