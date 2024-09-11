
# include("../plots_default.jl")
using PyPlot
using FLOWFoil
using FLOWMath

include("../test/data/bodyofrevolutioncoords.jl")

#Plot
figure(; figsize=(6, 4))
ylim([-0.25, 1.5])
xlabel("x")
ylabel(L"V_s/V_\infty")

plot(x, r, "C2"; linewidth=2, label="geometry")

# p = plot(; xlabel="x", ylabel=L"V_s/V_\infty", ylimits=[-0.25, 1.5])
# plot!(x, r; linecolor=mycolors[3], label="geometry")

# plot!(xsmooth, rsmooth; linestyle=:dash, label="smoothed geometry")

# mesh = [FLOWFoil.generate_axisym_mesh(xsmooth, rsmooth; bodyofrevolution=true)]
mesh = [FLOWFoil.generate_axisym_mesh(x, r; bodyofrevolution=true)]

problem = FLOWFoil.Problem(mesh; axisymmetric=true, viscous=false)

solution = FLOWFoil.solve(problem)

# get surface velocity at control points
cpx = [mesh[1].panels[i].controlpoint[1] for i in 1:length(solution.panelgammas)]
surface_velocity = solution.panelgammas

# experimental data
# plot!(
#    Vs_over_Vinf_x,
#    Vs_over_Vinf_vs;
#    seriestype=:scatter,
#    markeralpha=0,
#    markerstrokealpha=1,
#    markerstrokecolor=mycolors[1],
#    label="Experimental",
#)

plot(Vs_over_Vinf_x, Vs_over_Vinf_vs, "oC1"; label="Experimental")
plot(cpx, surface_velocity, "C0"; linewidth=2, label="Panel Method")

legend()

#panel method solution
# plot!(
#     cpx,
#     surface_velocity;
#     linecolor=mycolors[1],
#     label="panel method",
#     background_color=:transparent,
# )

# savefig("surface_velocity.pdf")
savefig("./docs/src/bodyofrevolution.pdf"; bbox_inches="tight")
# savetightplot(p, "surface_velocity.pdf")
