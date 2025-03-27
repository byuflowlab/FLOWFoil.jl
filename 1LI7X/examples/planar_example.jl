using FLOWFoil

#= NOTE:
    Very few combinations of singlarities and orders and boundary conditions are implemented right now.  Currently, this example provides the only implemented combination.
=#

#---------------------------------#
#           Single Body           #
#---------------------------------#

### --- Define Coordinates --- ###
# SEE: src/airfoils/parameterizations/naca.jl
x, z = naca4() #use one of the airfoil parameter convenience functions, default is 2412
coordinates = [x z] #coordinates must be in this format.

### --- Define Freestream Parameters --- ###
flow_angle = [0.0; 5.0; 10.0] #angles of attack in degrees (must be a vector or float)
reynolds = [-1.0] #!nothing implemented for reynolds number yet
mach = [-1.0] #!nothing implemented for mach number yet

### --- Set Up Problem --- ###
# SEE: src/dispatch_types.jl

# - Choose Order of Singularity - #
singularity_order = Linear()

# - Choose Singularity Type - #
singularity_type = Vortex(singularity_order)

# - Choose Boundary Condition Type - #
boundary_condition = Dirichlet()

# - Define Solution Method (Problem Type) - #
method = PlanarProblem(singularity_type, boundary_condition)

# - Generate Problem Object - #
# SEE: src/problem.jl
problem = define_problem(method, coordinates, flow_angle, reynolds, mach)

### --- Set Up the Rest --- ###

# - Generate Panel Geometry - #
# SEE: src/panel.jl
panels = generate_panels(method, coordinates)

# - Generate Influence Mesh - #
# SEE: src/mesh.jl, and src/geometry.jl
mesh, TEmesh = generate_mesh(method, panels)

# - Assemble Linear System - #
# SEE: src/system.jl, and src/singularity.jl
system = generate_inviscid_system(method, panels, mesh, TEmesh)

### --- Solve and Post Process --- ###
# - Solve Linear System - #
# SEE: src/solve.jl
solution = solve(system)

# - Post Process to get coefficients and surface distributions - #
# SEE: src/post_process.jl
polar = post_process(method, problem, panels, mesh, solution)

cl = polar.lift #total lift
cd = polar.drag #total drag
cm = polar.moment #total moment
surface_velocity = polar.surface_velocity #indexed as [body, panel, AoA]
surface_pressure = polar.surface_pressure #indexed as [body, panel, AoA]

#---------------------------------#
#    Using Convenience Function   #
#---------------------------------#
# SEE: src/convenience_functions.jl

polar2 = solve(coordinates, flow_angle, method)
cl2 = polar2.lift #total lift
cd2 = polar2.drag #total drag
cm2 = polar2.moment #total moment
surface_velocity2 = polar2.surface_velocity #indexed as [body, panel, AoA]
surface_pressure2 = polar2.surface_pressure #indexed as [body, panel, AoA]
