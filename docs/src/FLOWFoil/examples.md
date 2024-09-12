# Advanced Examples: Verification and Validation

<!-- ## Mfoil: Single inviscid airfoil comparision to analytic solution -->

<!-- TODO: add comparison with Joukowsky airfoil used in XFoil paper: -->
<!-- ```julia -->
<!-- using FLOWFoil -->

<!-- center = [-0.1; 0.1] -->
<!-- R = 1.0 -->
<!-- alpha = 4.0 -->

<!-- # get joukowsky coordinates from AirfoilTools -->

<!-- # run joukowsky solution from AirfoilTools -->

<!-- # run analyze function -->

<!-- # plot and save comparision, but hide the code -->
<!-- ``` -->

<!-- TODO: add comparision figure for single joukowsky airfoil here -->

<!-- ## Mfoil: Multiple inviscid airfoil comparison to analytic solution -->

<!-- For a multi-element airfoil system, the procedure is identical to a single body system, except we input a vector of matrices for the coordinates of the various bodies. -->
<!-- For this case, we'll use data that comes from ["An Exact Test Case for the Plane Potential Flow About Two Adjacent Lifting Aerofoils" by B. R. Williams.](https://reports.aerade.cranfield.ac.uk/handle/1826.2/2993) -->

<!-- ```@example -->
<!-- using FLOWFoil -->

<!-- # SET UP GEOMETRY -->
<!-- af_geom_path = normpath(joinpath(splitdir(pathof(FLOWFoil))[1], "..", "docs", "src", "assets", "two_inviscid_airfoils.jl")) -->
<!-- include(af_geom_path) -->

<!-- alpha = 0.0 -->

<!-- outputs = analyze([[ximain etamain], [xiflap etaflap]], alpha; method=Mfoil(inviscid=true) -->

<!-- plot and save comparisons, hiding code -->

<!-- nothing #hide -->
<!-- ``` -->

<!-- We see excellent agreement with the analytical solution. -->
<!-- TODO: add generated plot here -->
<!-- ![]() -->


<!-- ## Advanced Viscous Airfoil Options -->

<!-- TODO: explain the various advanced options (ncrit, forced transition, etc) -->

<!-- --- -->

## Axisymmetric Body of Revolution

For this example, we use data from chapter 4 of ["Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems" by  R. I. Lewis](https://doi.org/10.1017/CBO9780511529542)

```@example axisym
using FLOWFoil

bor_path = normpath(joinpath(splitdir(pathof(FLOWFoil))[1], "..", "test", "data", "bodyofrevolutioncoords.jl"))
include(bor_path)

outputs = analyze(center_body_coordinates, 0.0; method = Lewis(body_of_revolution=[true]))

# plot # hide
include("../assets/plots_default.jl")
plot(xlabel=L"\frac{x}{c}", ylabel=L"\frac{V_s}{V_\infty}")
plot!(Vs_over_Vinf_x, Vs_over_Vinf_vs, seriestype=:scatter, label="Experimental Data",markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(0.5*(center_body_coordinates[1:end-1,1].+center_body_coordinates[2:end,1]), outputs.tangential_velocities[1], label="FLOWFoil") #hide
savefig("../assets/body_of_revolution.png") #hide
nothing #hide
```

In order to have the correct flags associated with the axisymmetric solver, we generate the mehs using the `generate_axisym_mesh` function.
In addition, since we are modeling a body of revolution, that is, we have an open geometry at the axis of rotation, we need to include the `bodyofrevolution` keyword argument.

![](assets/body_of_revolution.png)

## Axisymmetric Annular Airfoil (Duct)

If we define an airfoil shape in an axisymmetric scheme, we model an annular airfoil, or in other words, a duct.  To do so, we follow a similar procedure to bodies of revolution with the exception that we set `body_of_revolution=false`.

```@example axisym
using FLOWFoil

duct_path = normpath(joinpath(splitdir(pathof(FLOWFoil))[1], "..", "test", "data", "naca_662-015.jl"))
include(duct_path)

outputs = analyze(duct_coordinates, 0.0; method = Lewis(body_of_revolution=[false]))

# plot # hide
include("../assets/plots_default.jl")
plot(xlabel=L"\frac{x}{c}", ylabel=L"c_p")
plot!(pressurexupper, pressureupper, seriestype=:scatter, markershape=:utriangle, label="Experimental Nacelle", color=1, yflip=true, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(pressurexlower, pressurelower, seriestype=:scatter, markershape=:dtriangle, label="Experimental Casing", color=1, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(0.5*(duct_coordinates[1:end-1,1].+duct_coordinates[2:end,1]), outputs.surface_pressures[1], label="FLOWFoil",color=2) #hide
savefig("../assets/annular_airfoil.png") #hide
nothing #hide
```

As above, we plot experimental results along with our calculated values.

![](assets/annular_airfoil.png)

## Axisymmetric Mutli-element Systems

As an example of an multi-element axisymmetric system (such as that used for a ducted rotor), we will simply combine the two previous cases.
Note that we include the coordinates for the various bodies as a vector of matrices, and in this case, we need to indicate in the `Lewis` method fields which of the bodies is a body of revolution (currently, the method only works if the annular airfoil comes first).

```@example axisym
using FLOWFoil

outputs = analyze([duct_coordinates, center_body_coordinates], 0.0; method = Lewis(body_of_revolution=[false, true]))


# plot v # hide
include("../assets/plots_default.jl")
plot(xlabel=L"\frac{x}{c}", ylabel=L"\frac{V_s}{V_\infty}")
plot!(Vs_over_Vinf_x, Vs_over_Vinf_vs, seriestype=:scatter, label="Experimental Center Body",markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(0.5*(center_body_coordinates[1:end-1,1].+center_body_coordinates[2:end,1]), outputs.tangential_velocities[2], label="FLOWFoil Center Body with Duct Effects") #hide
savefig("../assets/duct_hub_vs.png") #hide

# plot cp # hide
include("../assets/plots_default.jl")
plot(xlabel=L"\frac{x}{c}", ylabel=L"c_p")
plot!(pressurexupper, pressureupper, seriestype=:scatter, markershape=:utriangle, label="Experimental Nacelle", color=1, yflip=true, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(pressurexlower, pressurelower, seriestype=:scatter, markershape=:dtriangle, label="Experimental Casing", color=1, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(0.5*(duct_coordinates[1:end-1,1].+duct_coordinates[2:end,1]), outputs.surface_pressures[1], label="FLOWFoil Duct with Center Body Effects",color=2) #hide
savefig("../assets/duct_hub_cp.png") #hide
nothing #hide
```

Plotting the geometry and the output velocities and pressures show expected behavior when combining these two cases.

![](assets/duct_hub_vs.png)

![](assets/duct_hub_cp.png)
