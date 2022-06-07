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
polar = FLOWFoil.inviscid_post(inviscid_solution, alpha)

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
