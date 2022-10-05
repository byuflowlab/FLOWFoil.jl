
using FLOWFoil
using PyPlot

# af Geometry
# use a simple naca airfoil and rotate and scale and transform

xl, xu, yl, yu = naca4(2.0, 4.0, 20.0; split=true)
af = generate_mesh([xu yu]; thinaf=true)

# Plot Geometry
figure(1; figsize=(12, 5))
clf()
axis("equal")
# af
plot(getindex.(af.nodes, 1), getindex.(af.nodes, 2), "-"; linewidth=3)

axis("off")
savefig("ductgeom.pdf"; bbox_inches="tight")
axis("on")

# Define Grid Ranges
xrange = [
    minimum(getindex.(af.nodes, 1)) - af.chord / 2.0
    maximum(getindex.(af.nodes, 2)) + af.chord
]
zrange = [-0.5; 0.5]
Nz = 50
#make sure grid is square (not sure if that's necessary)
Nx = round(Int, Nz * (xrange[2] - xrange[1]) / (zrange[2] - zrange[1]))
xgp = range(xrange[1], xrange[2]; length=Nx)
zgp = range(zrange[1], zrange[2]; length=Nz)

# Grid Points (note, this grabs points outside the grid area, since only rectangles can be defined)
# plot(repeat(xgp; inner=(1, Nz)), repeat(zgp; inner=(1, Nx))', ".C2";)
#####################################
# Solve System
problem = Problem([af]; viscous=false)
solution = solve(problem)
# Solve stream values on grid
xg, zg, sg = calculate_stream_grid(problem, solution, xrange, zrange; Nx=Nx, Nz=Nz)

# Plot stream contours
contour(xg, zg, sg, 50; linestyles="-", colors="C2")

# plot([0.9; 1.1], [0.0; 0.0], "--C3"; linewidth=3)

# for slides
axis("off")
savefig("ductstream.pdf"; bbox_inches="tight")

legend()

