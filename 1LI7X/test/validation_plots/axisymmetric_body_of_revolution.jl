using FLOWFoil
include("../../plots_default.jl")
pyplot()

include("../data/bodyofrevolutioncoords.jl")

#Plot
plot(; xlabel=L"x", ylabel=L"V_s/V_\infty")

# plot!(x, r, "C2",linewidth=2, label="Geometry")
method = AxisymmetricProblem(Vortex(Constant()), Neumann(), [true])
flow_angle = [0.0]
reynolds = [-1.0]
mach = [-1.0]
coordinates = [x r]

# Generate Problem Object
problem = define_problem(method, coordinates, flow_angle, reynolds, mach)

# Generate Panel Geometry
panels = generate_panels(method, coordinates)

# Generate Influence Mesh
mesh = generate_mesh(method, panels)

# Assemble Linear System
system = generate_inviscid_system(method, panels, mesh)

# Solve Linear System
solution = solve(system)

# Post Process Solution
polar = post_process(method, problem, panels, mesh, solution)

# experimental data
plot!(
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    seriestype=:scatter,
    markershape=:utriangle,
    markercolor=mycolors[2],
    label="Experimental",
)

# panel method solution
plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_velocity[1, :, 1];
    linecolor=mycolors[1],
    label="FLOWFoil",
)

savefig("axisymmetric_body_of_revolution_validation_plot.pdf")
