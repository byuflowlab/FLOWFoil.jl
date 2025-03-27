using FLOWFoil
using Xfoil
include("../../plots_default.jl")
pyplot()

x, z = FLOWFoil.naca4(; blunt_TE=true)
# x, z = FLOWFoil.naca4()

N = length(x)

coordinates = ([x z], [x .+ 1000.0 z])

alpha = 0.0

Xfoil.set_coordinates(reverse(x), reverse(z))
cl, cm = Xfoil.solve_alpha(alpha)
g = Xfoil.get_globals()

xxf = g.x[1:N]
cpxf = g.cpi[1:N]
vxf = g.qinv[1:N]

polar = FLOWFoil.solve(coordinates, [alpha], PlanarProblem(Vortex(Linear()), Dirichlet()))

# PLOT
plot(; yflip=true, xlabel=L"\frac{x}{c}", ylabel=L"c_p")

plot!(xxf[1:end], cpxf[1:end]; label="Xfoil")

plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_pressure[1, :, 1];
    linestyle=:dash,
    linewidth=2,
    label="FLOWFoil",
)

# plot!(
#     polar.xsmooth[2, :, 1],
#     polar.surface_pressure[2, :, 1];
#     linestyle=:dash,
#     linewidth=2,
#     label="FLOWFoil",
# )

savefig("xfoil_comparison_plot_cp.pdf")

plot(; yflip=true, xlabel=L"\frac{x}{c}", ylabel=L"c_p")

plot!(xxf[1:end], vxf[1:end]; label="Xfoil")

plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_velocity[1, :, 1];
    linestyle=:dash,
    linewidth=2,
    label="FLOWFoil",
)

savefig("xfoil_comparison_plot_vt.pdf")
