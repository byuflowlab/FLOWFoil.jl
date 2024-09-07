using FLOWFoil
using PyPlot

figure(3; figsize=(6, 4))
clf()

include("../test/data/naca_662-015.jl")
duct = FLOWFoil.generate_axisym_mesh(x, r; bodyofrevolution=false)
plot(x, r; label="Annular Airfoil")

include("../test/data/bodyofrevolutioncoords.jl")
hub = FLOWFoil.generate_axisym_mesh(x, r; bodyofrevolution=true)
plot(x, r; label="Body of Revolution")
axis("equal")
legend()
savefig("ducthubgeom.pdf"; bbox_inches="tight")

meshes = [duct; hub]
problem = FLOWFoil.Problem([duct; hub]; axisymmetric=true, viscous=false)

solution = FLOWFoil.solve(problem)

# Post Processing

# get surface velocity at control points
# cpx = [mesh[1].panels[i].controlpoint[1] for i in 1:length(solution.panelgammas)]
cpx = [(p -> p.controlpoint[1]).(duct.panels); (p -> p.controlpoint[1]).(hub.panels)]

#surface velocities
gammas = solution.panelgammas

# surface_velocity = FLOWFoil.axisymmetric_surface_pressure(solution)
cp = 1.0 .- gammas .^ 2

using PyPlot
figure(1; figsize=(6, 4))
clf()
ylim([1.0, -1.2])
xlabel("x")
ylabel(L"c_P")

plot(pressurexlower, pressurelower, "oC2"; label="Experimental Inner Surface")
plot(pressurexupper, pressureupper, "^C2"; label="Experimental Outer Surface")
plot(cpx[1:40], cp[1:40]; label="Annular Airfoil")
plot(cpx[41:end], cp[41:end]; label="Body of Revolution")
legend()

savefig("multi_body_cp.pdf"; bbox_inches="tight")

figure(2; figsize=(6, 4))
clf()
ylim([-1.5, 1.5])
xlabel("x")
ylabel(L"V_s/V_\infty")
plot(Vs_over_Vinf_x, Vs_over_Vinf_vs, "oC2"; label="Experimental Body of Revolution")
plot(cpx[1:40], gammas[1:40]; label="Annular Airfoil")
plot(cpx[41:end], gammas[41:end]; label="Body of Revolution")
legend()

savefig("multi_body_vs.pdf"; bbox_inches="tight")

gamma_duct = FLOWFoil.get_mesh_gammas(gammas, meshes, 1)
gamma_hub = FLOWFoil.get_mesh_gammas(gammas, meshes, 2)
