using FLOWFoil

include("plots_default.jl")
pyplot()

# ANALYTIC SOLUTION
center = [-0.1; 0.1]
R = 1.0
alpha = 4.0
U = 1.0

xj, yj, vj, cpj, clj = FLOWFoil.joukowskysurface(center, R, alpha, U; N=120)

plot(xj, vj; label="Joukowsky")

# FLOWFOIL SOLUTION

meshes = [FLOWFoil.generate_mesh([xj yj])]
problem = FLOWFoil.Problem(meshes, alpha; viscous=false)

inviscid_solution = FLOWFoil.solve(problem)

polar = FLOWFoil.get_planar_polar(inviscid_solution, alpha)

plot!(xj, polar.surface_velocity; linestyle=:dash, linewidth=2, label="FLOWFoil")

savefig("single_airfoil_validation.jpg")
