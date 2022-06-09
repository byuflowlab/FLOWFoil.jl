# Quick Start

FLOWFoil is structured as follows:
 - The user generates a `Problem`.
 - The `Problem` is solved, generating a `Solution`.
 - The `Solution` is post-processed, generating a `Polar`.

There are also some included convenience functions for geometry generation and manipulation.

## Single Airfoil Inviscid Solution

Let's first look at the simplest case, a single inviscid airfoil.
We'll first set up the geometry, then defin the problem, then solve the problem, then post process it, and finally plot some of the outputs.

```@example
using FLOWFoil
using PyPlot

## -- SET UP GEOMETRY
# arbitrarily pick some joukowsky airfoil parameters
center = [-0.1; 0.1]
radius = 1.0
num_nodes = 160

# set freestream to unity
vinf = 1.0
re = 1.0

# arbitrarily pick an angle of attack
alpha = 4.0

# get airfoil coordinates for joukowsky airfoil
x, z = FLOWFoil.joukowsky(center, radius; N=num_nodes)

# get analytic joukowsky solution for later
vj, cpj, clj = FLOWFoil.joukowskyflow(center, radius, alpha, vinf; N=num_nodes)

# generate mesh object
meshes = [FLOWFoil.generate_mesh([x z])]

## -- DEFINE PROBLEM
problem = FLOWFoil.Problem(meshes, alpha, re; viscous=false)

## -- SOLVE PROBLEM
inviscid_solution = FLOWFoil.solve(problem)

## -- POST PROCESS SOLUTION
polar = FLOWFoil.inviscid_polar(inviscid_solution, alpha)

## -- PLOT
figure(1; figsize=(12, 3))

subplot(131)
plot(x, z)
axis("equal")
axis("off")

subplot(132)
xlabel("x")
ylabel(L"\frac{V_T}{V_\infty}")
plot(x, vj; label="Joukowsky")
plot(x, polar.surfacevelocity, "--"; linewidth=2, label="FLOWFoil")
legend()

subplot(133)
xlabel("x")
ylabel(L"c_p")
plot(x, cpj; label="Joukowsky")
plot(x, polar.surfacepressure, "--"; linewidth=2, label="FLOWFoil")
ylim(1.0, -1.75)
legend()
```

![](joukowsky.jpg)



## Multiple Airfoil Inviscid Solution

For a multi-element airfoil system, the procedure is identical, except an array of meshes is used for the problem definition.
For this case, we'll use data that comes from "An Exact Test Case for the Plane Potential Flow About Two Adjacent Lifting Aerofoils" by B. R. Williams.

```@example
using FLOWFoil
using PyPlot

## -- SET UP GEOMETRY
include("two_inviscid_airfoils.jl")

# set freestream to unity
vinf = 1.0
re = 1.0

# arbitrarily pick an angle of attack
alpha = 0.0

# generate mesh object
meshes = [FLOWFoil.generate_mesh([ximain etamain]); FLOWFoil.generate_mesh([xiflap etaflap])]

## -- DEFINE PROBLEM
problem = FLOWFoil.Problem(meshes, alpha, re; viscous=false)

## -- SOLVE PROBLEM
inviscid_solution = FLOWFoil.solve(problem)

## -- POST PROCESS SOLUTION
polar = FLOWFoil.inviscid_polar(inviscid_solution, alpha)

## -- PLOT
figure(2; figsize=(9, 3))
clf()
subplot(121)
plot(ximain, etamain, label="Main Airfoil")
plot(xiflap, etaflap, label="Flap Airfoil")
axis("equal")
axis("off")

subplot(122)
xlabel("x")
ylabel(L"c_p")
plot(ximain, cpmain; label="Analytic")
plot(xiflap, cpflap)
plot(ximain, polar.surfacepressure[1:length(ximain)], ".C0"; linewidth=2, label="FLOWFoil")
plot(xiflap, polar.surfacepressure[length(ximain)+1:end], ".C1"; linewidth=2)
ylim(maximum(polar.surfacepressure)*1.1, minimum(polar.surfacepressure)*1.1)
legend()
```

![](two_inviscid_airfoils.jpg)
