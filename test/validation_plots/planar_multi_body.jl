using FLOWFoil
include("../../plots_default.jl")
pyplot()

## -- SET UP GEOMETRY
include("../data/planar_multi_body_data.jl")

coordinates = ([ximain etamain], [xiflap etaflap])

polar = FLOWFoil.solve(coordinates, PlanarProblem(Vortex(Linear()), Dirichlet()))

## -- PLOT
plot(; yflip=true, xlabel=L"x", ylabel=L"c_p")
plot!(
    ximain,
    cpmain;
    seriestype=:scatter,
    markershape=:utriangle,
    markercolor=mycolors[1],
    label="Analytic Main Airfoil",
)
plot!(
    xiflap,
    cpflap;
    seriestype=:scatter,
    markershape=:utriangle,
    markercolor=mycolors[2],
    label="Analytic Flap Airfoil",
)

plot!(
    polar.xsmooth[1, :, 1],
    polar.surface_pressure[1, :, 1];
    color=mycolors[2],
    label="FLOWFoil Main Airfoil",
)

plot!(
    polar.xsmooth[2, :, 1],
    polar.surface_pressure[2, :, 1];
    color=mycolors[1],
    label="FLOWFoil Flap Airfoil",
)

savefig("planar_multi_body_validation_plot.pdf")
