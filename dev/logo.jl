using FLOWFoil
using PyPlot

figure(1; figsize=(4, 4))
clf()

# x, y = naca4(2.0, 5.0, 24.0)
# m1 = FLOWFoil.generate_mesh(x, y)

# x, y = naca4(6.0, 4.0, 24.0)
# m2 = FLOWFoil.generate_mesh(x, y)

# x, y = naca4(9.0, 5.0, 24.0)
# m3 = FLOWFoil.generate_mesh(x, y)

# inviscid_problem = FLOWFoil.Problem([m1; m3], 4.0, 0.0; viscous=false)
# inviscid_solution = FLOWFoil.solve(inviscid_problem)
# polar = FLOWFoil.inviscid_polar(inviscid_solution, 0.0)
# x, z, stream = FLOWFoil.calculate_stream_grid(
#     inviscid_problem, inviscid_solution, [-1.0; 1.9], [-1.15; 0.65]
# )

# contour(x, z, stream, 50; linestyles="-", colors=[(0.251, 0.388, 0.847)])

xu, xl, yu, yl = naca4(12.0, 5.0, 44.0; split=true)
fill_between(xu, yu, reverse(yl); color=(0.22, 0.596, 0.149))

xu, xl, yu, yl = naca4(12.0, 6.0, 44.0; split=true)
xu .+= 0.75
xl .+= 0.75
yu .-= 0.75
yl .-= 0.75
fill_between(xu, yu, reverse(yl); color=(0.584, 0.345, 0.698))

xu, xl, yu, yl = naca4(8.0, 5.0, 44.0; split=true)
xu .-= 0.9
xl .-= 0.9
yu .-= 0.75
yl .-= 0.75
fill_between(xu, yu, reverse(yl); color=(0.796, 0.235, 0.2))

axis("equal")
axis("off")

savefig("src/assets/logo.png"; bbox_inches="tight")
