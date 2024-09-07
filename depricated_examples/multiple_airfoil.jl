using FLOWFoil

## -- SET UP GEOMETRY
include("../docs/src/two_inviscid_airfoils.jl")

# arbitrarily pick an angle of attack
alpha = 0.0

# generate mesh object
meshes = [
    FLOWFoil.generate_mesh([ximain etamain])
    FLOWFoil.generate_mesh([xiflap etaflap])
]

## -- DEFINE PROBLEM
problem = FLOWFoil.Problem(meshes; viscous=false)

## -- SOLVE PROBLEM
inviscid_solution = FLOWFoil.solve(problem)

## -- POST PROCESS SOLUTION
polar = FLOWFoil.get_planar_polar(inviscid_solution, alpha)

## -- PLOT
using PyPlot
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
plot(xiflap, cpflap, "C0")
plot(ximain, polar.surface_pressure[1:length(ximain)], "--C1"; linewidth=2, label="FLOWFoil")
plot(xiflap, polar.surface_pressure[length(ximain)+1:end], "--C1"; linewidth=2)
ylim(maximum(polar.surface_pressure)*1.1, minimum(polar.surface_pressure)*1.1)
legend()

# include("../plots_default.jl")
# af = plot(; aspectratio=:equal, axis=false)
# plot!(af, ximain, etamain; label="Main Airfoil")
# plot!(af, xiflap, etaflap; label="Flap Airfoil")

# savefig(af, "./../../Writing/whisper_sttr/two_airfoils.tikz")

# cps = plot(; yflip=true, xlabel="x-coordinate", ylabel="Surface Pressure Coefficient")
# plot!(cps, ximain, cpmain; color=mycolors[1], label="Analytic")
# plot!(cps, xiflap, cpflap; color=mycolors[1], label="")
# plot!(
#     cps,
#     ximain,
#     polar.surface_pressure[1:length(ximain)];
#     # seriestype=:scatter,
#     linestyle=:dash,
#     linewidth=2,
#     color=mycolors[2],
#     label="FLOWFoil",
# )
# plot!(
#     cps,
#     xiflap,
#     polar.surface_pressure[(length(ximain) + 1):end];
#     # seriestype=:scatter,
#     linestyle=:dash,
#     linewidth=2,
#     color=mycolors[2],
#     label="",
# )

# savefig(cps, "./../../Writing/whisper_sttr/multi_airfoil_validation.tikz")
