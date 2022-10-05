using FLOWFoil
using PyPlot
using Xfoil

# ANALYTIC SOLUTION
center = [-0.1; 0.1]
R = 1.0
alpha = 4.0
U = 1.0

xj, yj, vj, cpj, clj = FLOWFoil.joukowskysurface(center, R, alpha, U; N=120)

figure(1; figsize=(8, 3))
clf()

subplot(121)
xlabel("x")
ylabel(L"\frac{V_T}{V_\infty}")
plot(xj, vj; label="Joukowsky")

subplot(122)
xlabel("x")
ylabel(L"c_p")
plot(xj, cpj; label="Joukowsky")
ylim(1.0, -1.75)

# FLOWFOIL SOLUTION
re = 1.0

meshes = [FLOWFoil.generate_mesh([xj yj])]
problem = FLOWFoil.Problem(meshes, alpha, re; viscous=false)

inviscid_solution = FLOWFoil.solve(problem)

polar = FLOWFoil.inviscid_post(inviscid_solution, alpha)

subplot(121)
plot(xj, polar.surfacevelocity, "--"; linewidth=2, label="FLOWFoil")
legend()

subplot(122)
plot(xj, polar.surfacepressure, "--"; linewidth=2, label="FLOWFoil")
legend()

savefig("joukowsky_comp.pdf"; bbox_inches="tight")
