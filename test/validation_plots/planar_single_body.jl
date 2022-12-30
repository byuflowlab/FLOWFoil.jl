using FLOWFoil
using Xfoil

include("../../plots_default.jl")
pyplot()

# ANALYTIC SOLUTION
center = [-0.1; 0.1]
R = 1.0
alpha = 4.0
U = 1.0

xj, yj = FLOWFoil.joukowsky(center, R; N=120)
vsurf, cpsurf, cl = FLOWFoil.joukowskyflow(center, R, alpha, U; N=120)

# FLOWFOIL SOLUTION

polar = FLOWFoil.solve([xj yj], [alpha], PlanarProblem(Vortex(Linear()), Dirichlet()))

# polar = FLOWFoil.solve(
#     [xj[1:(end - 1)] yj[1:(end - 1)]], [alpha], PlanarProblem(Vortex(Linear()), Dirichlet())
# )

# PLOT
plot(; yflip=true, xlabel=L"\frac{x}{c}", ylabel=L"c_p")

# plot!(xj[2:end], cpsurf[2:end]; label="Joukowsky")
plot!(xj[1:end], cpsurf[1:end]; label="Joukowsky")

plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_pressure[1, :, 1];
    linestyle=:dash,
    linewidth=2,
    label="FLOWFoil",
)

savefig("planar_single_body_validation_plot.pdf")
