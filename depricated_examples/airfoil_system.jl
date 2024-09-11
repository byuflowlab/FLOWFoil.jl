using FLOWFoil
using PyPlot

# ANALYTIC SOLUTION
center = [-0.1; 0.1]
R = 1.0
alpha = 0.0
U = 1.0

xj, yj, vj, cpj, clj = FLOWFoil.joukowskysurface(center, R, alpha, U; N=120)
xj .-= minimum(xj)
xscale = maximum(xj) - minimum(xj)
xj /= xscale
yj /= xscale

# xj = [1.0; 0.5; 0.0; 0.5; 1.0]
# yj = [0.0; -0.5; 0.0; 0.5; 0.0]

# FLOWFOIL SOLUTION
re = 1.0
scales = [1.0; 1.0]
angles = [0.0; 0.0]
locations = [[0.0, 0.0], [2.0, 0.0]]
meshes = [FLOWFoil.generate_mesh([xj yj]); FLOWFoil.generate_mesh([xj yj])]
FLOWFoil.position_coordinates!(meshes, angles, scales, locations)

figure(1; figsize=(9, 3))
clf()

subplot(131)
for i in 1:length(meshes)
    plot(getindex.(meshes[i].airfoil_nodes, 1), getindex.(meshes[i].airfoil_nodes, 2))
end
axis("equal")

problem = FLOWFoil.Problem(meshes, alpha, re; viscous=false)

inviscid_solution = FLOWFoil.solve(problem)

polar = FLOWFoil.inviscid_post(inviscid_solution, alpha)

subplot(132)
xlabel("x")
ylabel(L"\frac{V_T}{V_\infty}")
offset = 0
for i in 1:length(meshes)
    plot(
        getindex.(meshes[i].airfoil_nodes, 1),
        polar.surfacevelocity[(1 + offset):(inviscid_solution.Ns[i] + offset)],
    )
    global offset += inviscid_solution.Ns[i]
end

subplot(133)
xlabel("x")
ylabel(L"c_p")
offset = 0
for i in 1:length(meshes)
    plot(
        getindex.(meshes[i].airfoil_nodes, 1),
        polar.surfacepressure[(1 + offset):(inviscid_solution.Ns[i] + offset)],
    )
    global offset += inviscid_solution.Ns[i]
end
ylim(1.0, -1.25)

savefig("twofoiltest.pdf"; bbox_inches="tight")
