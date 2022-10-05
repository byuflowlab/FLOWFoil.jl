using FLOWFoil
using PyPlot

# Wall Geometry
# use a simple naca airfoil and rotate and scale and transform

xwall, ywall = naca4(2.0, 4.0, 30.0)

scale = 1.0
location = [0.0, 0.75]
angle = 0.0

xwall, ywall = position_coordinates([xwall -ywall], scale, angle, location)

wall = generate_mesh([xwall ywall])
# wall_mirror = generate_mesh([xwall -ywall])

# Hub Geometry
# use a symmetric naca airfoil for now, worry about super blunt TE later
xl, xu, yl, yu = naca4(0.0, 0.0, 40.0; split=true)
hub = generate_mesh([xu yu]; axisymmetric=true)
# x, y = naca4(0.0, 0.0, 40.0)
# hub = generate_mesh([x y])



include("../../DuctTAPE/data/dfdc/dstestr2_case.jl")

wall = generate_mesh([ductx ductr])
x = [reverse(hubx); hubx[2:end]]
y = [reverse(-hubr); hubr[2:end]]
x = hubx
y = hubr
hub = generate_mesh([x y])

# Plot Geometry
figure(1; figsize=(12, 5))
clf()
axis("equal")
# Wall
plot(getindex.(wall.nodes, 1), getindex.(wall.nodes, 2); linewidth=3)
# plot(getindex.(wall_mirror.nodes, 1), getindex.(wall_mirror.nodes, 2), "C0"; linewidth=3)
# Hub
plot(getindex.(hub.nodes, 1), getindex.(hub.nodes, 2); linewidth=3)

axis("off")
savefig("ductgeom.pdf",bbox_inches="tight")
axis("on")

# Define Grid Ranges
xrange = [
    minimum(getindex.(wall.nodes, 1)) - wall.chord / 2.0
    maximum(getindex.(wall.nodes, 2)) + wall.chord
]
zrange = [-0.; wall.nodes[round(Int,length(wall.nodes)/2)][2]]
Nz = 25
#make sure grid is square (not sure if that's necessary)
Nx = round(Int, Nz * (xrange[2] - xrange[1]) / (zrange[2] - zrange[1]))
xgp = range(xrange[1], xrange[2]; length=Nx)
zgp = range(zrange[1], zrange[2]; length=Nz)

# Plot Grid
# Grid Boundaries
# plot(
#     [xrange[1]; xrange[1]],
#     [zrange[1]; zrange[2]],
#     "--k";
#     linewidth=4,
#     label="Solver Bounds",
# )
# plot([xrange[2]; xrange[2]], [zrange[1]; wall.nodes[1][2]], "--k"; linewidth=4)
# plot(
#     [wall.nodes[1][1]; xrange[2]], [wall.nodes[1][2]; wall.nodes[1][2]], "--k"; linewidth=4
# )
# plot([xrange[1]; location[1][1]], [zrange[2]; zrange[2]], "--k"; linewidth=4)
# plot([xrange[1]; 0.0], [zrange[1]; zrange[1]], "--k"; linewidth=4)
# plot([hub.nodes[end][1]; xrange[2]], [zrange[1]; zrange[1]], "--k"; linewidth=4)

# Grid Points (note, this grabs points outside the grid area, since only rectangles can be defined)
plot(repeat(xgp; inner=(1, Nz)), repeat(zgp; inner=(1, Nx))', ".C2";)
#####################################
# Solve System
problem = Problem([wall, hub]; viscous=false)
solution = solve(problem)
# Solve stream values on grid
xg, zg, sg = calculate_stream_grid(problem, solution, xrange, zrange; Nx=Nx, Nz=Nz)

# Plot stream contours
contour(xg, zg, sg, 50; linestyles="-", colors="C2")

# for slides
axis("off")
savefig("ductstream.pdf"; bbox_inches="tight")

legend()
