
using FLOWFoil
using PyPlot

# SET UP GEOMETRY
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

# DEFINE PROBLEM
problem = FLOWFoil.Problem(meshes, alpha, re; viscous=false)

# SOLVE PROBLEM
inviscid_solution = FLOWFoil.solve(problem)

# POST PROCESS SOLUTION
polar = FLOWFoil.inviscid_polar(inviscid_solution, alpha)

# PLOT
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

savefig("docs/src/joukowsky.jpg", bbox_inches="tight")



########
#MULTI AIRFOIL EXAMPLE
########

using FLOWFoil
using PyPlot

## -- SET UP GEOMETRY
include("../docs/src/two_inviscid_airfoils.jl")

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
legend()

subplot(122)
xlabel("x")
ylabel(L"c_p")
plot(ximain, cpmain; label="Analytic")
plot(xiflap, cpflap)
plot(ximain, polar.surfacepressure[1:length(ximain)], ".C0"; linewidth=2, label="FLOWFoil")
plot(xiflap, polar.surfacepressure[length(ximain)+1:end], ".C1"; linewidth=2)
ylim(maximum(polar.surfacepressure)*1.1, minimum(polar.surfacepressure)*1.1)
legend()

savefig("docs/src/two_inviscid_airfoils.jpg",bbox_inches="tight")
