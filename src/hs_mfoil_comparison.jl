using Plots
using FLOWFoil

# ########## Validate with Joukowsky #################
import FLOWFoil.AirfoilTools as at


include("C:\\Users\\nlehn\\my_defaults\\my_plots_default.jl")
# default(; size=subfig_size(; n_cols=2), markerstrokealpha=0, markersize=2)

# - Parameters - #
center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0
Vinf = 1.0

# - Joukowsky Geometry - #
x, y = at.joukowsky(center, radius)
xy = [x y]

# - Surface Values - #
surface_velocity, surface_pressure_coefficient, cl = at.joukowsky_flow(
    center, radius, alpha, Vinf
)

# - Your Stuff - #

# - Plot Stuff - #
pl = plot(; xlabel="X", ylabel="CP", yflip=true)
plot!(
    pl,
    x[7:360],
    surface_pressure_coefficient[7:360];
    linestyle=:dash,
    linewidth=2,
    label="Analytic Solution",
    title = "Coefficient of Pressure"
)

mf_results = FLOWFoil.analyze(xy, [alpha], [1e-6], [0.0]; method=Mfoil())
hs_results = FLOWFoil.analyze(xy, [alpha], [1e-6], [0.0]; method=HessSmith())

plot!(pl, x[10:350], mf_results.surface_pressures[1][10:350], label="Mfoil")
plot!(pl, x[10:350], hs_results.surface_pressures[1][10:350], label="Hess Smith")

display(pl)
# savefig(pl, "C:\\Users\\nlehn\\HSPanel\\2D\\HSvsMFoil_comparison.png")
