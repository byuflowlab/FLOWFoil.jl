using FLOWFoil
include("../../plots_default.jl")
pyplot()

include("../data/bodyofrevolutioncoords.jl")
xhub = x
rhub = r

include("../data/naca_662-015.jl")
xduct = x
rduct = r

method = AxisymmetricProblem(Vortex(Constant()), Neumann(), [true, false])
flow_angle = [0.0]
reynolds = [-1.0]
mach = [-1.0]
coordinates = ([xhub rhub], [xduct rduct])

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
plot(
    pressurexlower,
    pressurelower;
    seriestype=:scatter,
    markershape=:dtriangle,
    markercolor=mycolors[3],
    label="Experimental Inner Surface",
)
plot!(
    pressurexupper,
    pressureupper;
    seriestype=:scatter,
    markershape=:utriangle,
    markercolor=mycolors[3],
    label="Experimental Outer Surface",
)

# panel method solution
plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_pressure[1, :, 1];
    linecolor=mycolors[1],
    label="FLOWFoil",
)

plot!(
    polar.xsmooth[2, :, 1],
    polar.surface_pressure[2, :, 1];
    linecolor=mycolors[2],
    label="FLOWFoil",
)

savefig("axisymmetric_presssure_plot.pdf")

#Plot
plot(; xlabel=L"x", ylabel=L"V_s/V_\infty")

# experimental data
plot!(
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    seriestype=:scatter,
    markershape=:utriangle,
    markercolor=mycolors[3],
    label="Experimental",
)

# panel method solution
plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_velocity[1, :, 1];
    linecolor=mycolors[1],
    label="FLOWFoil",
)

# panel method solution
plot!(
    polar.xsmooth[2, :, 1],
    polar.surface_velocity[2, :, 1];
    linecolor=mycolors[2],
    label="FLOWFoil",
)

savefig("axisymmetric_velocity_plot.pdf")
