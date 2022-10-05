
include("../plots_default.jl")
using FLOWFoil
using FLOWMath

x = [1.0; 0.5; 0.0; 0.5; 1.0]
r = [0.0; -0.5; 0.0; 0.5; 0.0] .+ 1.5

x, r = naca4(N=100)
r .+= 0.5
p = plot(; xlabel="x/c", ylabel=L"c_P", yflip=true )#, aspectratio=:equal)

mesh = [FLOWFoil.generate_axisym_mesh(x, r; bodyofrevolution=false)]

problem = FLOWFoil.Problem(mesh; axisymmetric=true, viscous=false)

solution = FLOWFoil.solve(problem)

# get surface velocity at control points
cpx = [mesh[1].panels[i].controlpoint[1] for i in 1:(length(solution.panelgammas) - 1)]
# surface_velocity = FLOWFoil.axisymmetric_surface_pressure(solution)
cp = 1.0 .- solution.panelgammas .^ 2

#panel method solution
plot!(cpx, cp; label="panel method")

savefig("test_axisym_af.pdf")
