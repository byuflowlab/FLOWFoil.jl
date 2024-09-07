
# include("../plots_default.jl")
using FLOWFoil
using FLOWMath

include("../test/data/naca_662-015.jl")

#Plot
# p = plot(; xlabel="x/c", ylabel=L"c_P", ylimit=(-1.25, 1.0), yflip=true)#, aspectratio=:equal)
# plot!(x, r; label="geometry")

mesh = [FLOWFoil.generate_axisym_mesh(x, r; bodyofrevolution=false)]

problem = FLOWFoil.Problem(mesh; axisymmetric=true, viscous=false)

solution = FLOWFoil.solve(problem)

gammas = solution.panelgammas
# get surface velocity at control points
cpx = [mesh[1].panels[i].controlpoint[1] for i in 1:length(gammas)]
# surface_velocity = FLOWFoil.axisymmetric_surface_pressure(solution)
cp = 1.0 .- gammas .^ 2

#experimental data
#plot!(
#    pressurexlower,
#    pressurelower;
#    seriestype=:scatter,
#    markeralpha=0,
#    markerstrokecolor=mycolors[1],
#    markerstrokealpha=1,
#    label="Experimental (inside)",
#)
#plot!(
#    pressurexupper,
#    pressureupper;
#    seriestype=:scatter,
#    markeralpha=0,
#    markerstrokecolor=mycolors[2],
#    markerstrokealpha=1,
#    label="Experimental (outside)",
#)

##panel method solution
#plot!(cpx, cp; label="panel method", background_color=:transparent)

# savefig("annular_airfoil.pdf")

using PyPlot
figure(1; figsize=(6, 4))
clf()
ylim([0.8, -0.8])
xlabel("x")
ylabel(L"c_P")

plot(pressurexlower, pressurelower, "oC1"; label="Experimental Inner Surface")
plot(pressurexupper, pressureupper, "^C1"; label="Experimental Outer Surface")
plot(cpx, cp, "C0"; label="Panel Method")
legend()

savefig("annular_airfoil.pdf"; bbox_inches="tight")
