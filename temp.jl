using Plots, FLOWFoil, .AirfoilTools, DelimitedFiles, FLOWMath

using Plots.PlotMeasures

include("src\\AirfoilTools\\parameterizations\\naca_65series.jl")
include("src/AirfoilTools/airfoil_utilities.jl")
include("src/AirfoilTools/geometry_manipulations/utils.jl")
alpha_fig_30 = [-0.8181092834383961, 2.2024718495542945, 5.249261562686591, 10.213842991582608, 12.183544604784126, 14.232210133632053, 16.210366127523763, 21.134028353879238, 25.14022118773481]
cl_fig_30 = [0.15120659864412866, 0.33164844571495, 0.45534871678353084, 0.5870299232236553, 0.6426259306423745, 0.6872650606866014, 0.7245573339110494, 0.7848286138514461, 0.7514211285541424]

alpha_martensen =  alpha_fig_30
cl_martensen = similar(alpha_martensen) .= 0.0
x,y = naca65c(1.5, 1.92*10^-2, 1.0, "3-018")
stagger_martensen = similar(alpha_martensen) .= 0.0
inflow_angle = 45.0
solidity = 1.0
plot(x,y, aspectratio=1)
#compute stagger angle from inflow and outflow
#=
for i = 1:length(alpha_martensen)
    stagger_martensen[i] = inflow_angle - alpha_martensen[i]
    method1 = FLOWFoil.Martensen(true, solidity, stagger_martensen[i], 10.0, 100.0, false)
    sol = FLOWFoil.analyze(x,y, [inflow_angle], method = method1)
    cl_martensen[i] = sol.lift_coefficients[1][1]
end
println(cl_martensen)
=#
#=
method1 = FLOWFoil.Martensen(true, 1.0, 45.0, 10.0, 100.0, false)
sol = FLOWFoil.analyze(x,y, alpha_martensen, method = method1)
cl_martensen = sol.lift_coefficients[1][:]


scatter(alpha_fig_30, cl_fig_30, label = "NACA Memorandum - NACA 65-410 σ=0.75")
plot!(alpha_martensen, cl_martensen, label = "Martensen's Method - NACA 65-65-410 σ=0.75", grid = false, xlabel=("α"), left_margin = 8mm)
annotate!(-3.0, 0.6, text("Cl", 10))
=#
