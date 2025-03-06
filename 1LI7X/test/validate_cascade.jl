using FLOWFoil
include("cascade_data.jl")

coordinates = [x y]

pt = PeriodicProblem(Vortex(Constant()), Neumann(), 0.900364, 0.0)

problem = define_problem(pt, coordinates, 35.0, -1.0, -1.0)

# Generate Panel Geometry
panels = generate_panels(pt, coordinates)

# Generate Influence Mesh
mesh = generate_mesh(pt, panels)

# Assemble Linear System
system = generate_inviscid_system(pt, panels, mesh)

solution = solve(system)

polar = post_process(pt, problem, panels, mesh, solution)
