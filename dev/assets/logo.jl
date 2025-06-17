using FLOWFoil
using PyPlot

figure(1; figsize=(4, 4))
clf()

x, y = naca4(12.0, 5.0, 44.0)
m1 = FLOWFoil.generate_mesh(x, y)

x, y = naca4(12.0, 6.0, 44.0)
m2 = FLOWFoil.generate_mesh(x .+ 0.75, y .- 0.75)

x, y = naca4(8.0, 5.0, 44.0)
m3 = FLOWFoil.generate_mesh(x .- 0.9, y .- 0.75)

inviscid_problem = FLOWFoil.Problem([m1; m3; m2], 4.0, 0.0; viscous=false)
inviscid_solution = FLOWFoil.solve(inviscid_problem)
polar = FLOWFoil.get_planar_polar(inviscid_solution, 0.0)
x, y, stream = FLOWFoil.calculate_stream_grid(
    inviscid_problem, inviscid_solution, [-1.0; 1.9], [-1.15; 0.65]
)

contour(x, y, stream, 35; linewidths=1.5, linestyles="-", colors=[(0.251, 0.388, 0.847)])

xu, xl, yu, yl = naca4(12.0, 5.0, 44.0; split=true)
fill_between(xu, yu, reverse(yl); color=(0.22, 0.596, 0.149), zorder=2)
plot(xu, yu; color=(0.22, 0.596, 0.149), linewidth=2, zorder=3)
plot(xl, yl; color=(0.22, 0.596, 0.149), linewidth=2, zorder=3)

xu, xl, yu, yl = naca4(12.0, 6.0, 44.0; split=true)
xu .+= 0.75
xl .+= 0.75
yu .-= 0.75
yl .-= 0.75
fill_between(xu, yu, reverse(yl); color=(0.584, 0.345, 0.698), zorder=2)
plot(xu, yu; color=(0.584, 0.345, 0.698), linewidth=2, zorder=3)
plot(xl, yl; color=(0.584, 0.345, 0.698), linewidth=2, zorder=3)

xu, xl, yu, yl = naca4(8.0, 5.0, 44.0; split=true)
xu .-= 0.9
xl .-= 0.9
yu .-= 0.75
yl .-= 0.75
fill_between(xu, yu, reverse(yl); color=(0.796, 0.235, 0.2), zorder=2)
plot(xu, yu; color=(0.796, 0.235, 0.2), linewidth=2, zorder=3)
plot(xl, yl; color=(0.796, 0.235, 0.2), linewidth=2, zorder=3)

axis("equal")
axis("off")

savefig("assets/logo.png"; bbox_inches="tight", transparent="true")
# savefig("logo.pdf"; bbox_inches="tight")
