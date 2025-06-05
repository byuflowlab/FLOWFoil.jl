# Advanced Examples: Verification and Validation

## Mfoil: Single inviscid airfoil comparision to analytic solution

TODO: add comparison with Joukowsky airfoil used in XFoil paper:
```@example Joukowsky
using FLOWFoil
using Plots
using .AirfoilTools

center = [-0.1; 0.1]
R = 1.0
alpha = 4.0
Vinf = 1.0 

# - Joukowsky Geometry - #
x, y = FLOWFoil.AirfoilTools.joukowsky(center, radius)

# - Surface Values - #
surface_velocity, surface_pressure_coefficient, cl = FLOWFoil.AirfoilTools.joukowsky_flow(
    center, radius, alpha, Vinf
)

# - Plot Stuff - #
pl = plot(; xlabel="x", ylabel="cp", yflip=true)
plot!(
    pl,
    x[7:340],
    surface_pressure_coefficient[7:340];
    linestyle=:dash,
    linewidth=2,
    label="Analytic Solution",
)

outputs = analyze([x,y], alpha; method="Mfoil")
plot!(pl, outputs.x[7:360], outputs.cp[7:360], label="Mfoil")
```

TODO: add comparision figure for single joukowsky airfoil here

## Mfoil: Multiple inviscid airfoil comparison to analytic solution

For a multi-element airfoil system, the procedure is identical to a single body system, except we input a vector of matrices for the coordinates of the various bodies.
For this case, we'll use data that comes from ["An Exact Test Case for the Plane Potential Flow About Two Adjacent Lifting Aerofoils" by B. R. Williams.](https://reports.aerade.cranfield.ac.uk/handle/1826.2/2993)

```@example Mfoil
using FLOWFoil

af_geom_path = normpath(joinpath(splitdir(pathof(FLOWFoil))[1], "..", "docs", "src", "assets", "two_inviscid_airfoils.jl"))
include(af_geom_path)

outputs = analyze([[ximain etamain], [xiflap etaflap]], [0.0]; method=Mfoil())

# TODO plot and save comparisons, hiding code

nothing #hide
```
![](assets/two_inviscid_airfoils.jpg)

We see excellent agreement with the analytical solution.

---

## Axisymmetric Body of Revolution

For this example, we use data from chapter 4 of ["Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems" by  R. I. Lewis](https://doi.org/10.1017/CBO9780511529542)

```@example axisym
using FLOWFoil

data_path = normpath(
    joinpath(
        splitdir(pathof(FLOWFoil))[1], "..", "test", "data", "bodyofrevolutioncoords.jl"
    ),
)
include(data_path)

outputs = analyze(center_body_coordinates, [0.0]; method=Lewis(; body_of_revolution=[true]))

nothing # hide

# plot # hide

# include("../assets/plots_default.jl") #hide
# plot(xlabel=L"\frac{x}{c}", ylabel=L"\frac{V_s}{V_\infty}") #hid
# plot!(Vs_over_Vinf_x, Vs_over_Vinf_vs, seriestype=:scatter, label="Experimental Data",markerstrokecolor=1, markercolor=1, markersize=4) #hide
# plot!(0.5*(center_body_coordinates[1:end-1,1].+center_body_coordinates[2:end,1]), outputs.vs, label="FLOWFoil") #hide
```
![](assets\bodyofrevolution.jpg)

## Axisymmetric Annular Airfoil (Duct)

If we define an airfoil shape in an axisymmetric scheme, we model an annular airfoil, or in other words, a duct.  To do so, we follow a similar procedure to bodies of revolution with the exception that we set `body_of_revolution=false`.

```@example axisym
using FLOWFoil

duct_path = normpath(
    joinpath(splitdir(pathof(FLOWFoil))[1], "..", "test", "data", "naca_662-015.jl")
)
include(duct_path)

outputs = analyze(duct_coordinates, [0.0]; method=Lewis(; body_of_revolution=[false]))

# plot # hide
plot(xlabel=L"\frac{x}{c}", ylabel=L"c_p") #hide
plot!(pressurexupper, pressureupper, seriestype=:scatter, markershape=:utriangle, label="Experimental Nacelle", color=1, yflip=true, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(pressurexlower, pressurelower, seriestype=:scatter, markershape=:dtriangle, label="Experimental Casing", color=1, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(0.5*(duct_coordinates[1:end-1,1].+duct_coordinates[2:end,1]), outputs.cp, label="FLOWFoil",color=2) #hide
```

As above, we plot experimental results along with our calculated values.


## Axisymmetric Mutli-element Systems

As an example of an multi-element axisymmetric system (such as that used for a ducted rotor), we will simply combine the two previous cases.
Note that we include the coordinates for the various bodies as a tuple of matrices, and in this case, we need to indicate in the `Lewis` method fields which of the bodies is a body of revolution (currently, the method only works if the annular airfoil comes first).

```@example axisym
using FLOWFoil

outputs = analyze(
    (duct_coordinates, center_body_coordinates), [0.0];
    method=Lewis(; body_of_revolution=[false, true]),
)

# plot v # hide
#=
plot(xlabel=L"\frac{x}{c}", ylabel=L"\frac{V_s}{V_\infty}") #hide
plot!(Vs_over_Vinf_x, Vs_over_Vinf_vs, seriestype=:scatter, label="Experimental Center Body",markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(0.5*(center_body_coordinates[1:end-1,1].+center_body_coordinates[2:end,1]), outputs.vs[2], label="FLOWFoil Center Body with Duct Effects") #hide
=#
```

```@example axisym
println(outputs)
# plot cp # hide
#=
plot(xlabel=L"\frac{x}{c}", ylabel=L"c_p") #hide
plot!(pressurexupper, pressureupper, seriestype=:scatter, markershape=:utriangle, label="Experimental Nacelle", color=1, yflip=true, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(pressurexlower, pressurelower, seriestype=:scatter, markershape=:dtriangle, label="Experimental Casing", color=1, markerstrokecolor=1, markercolor=1, markersize=4) #hide
plot!(0.5*(duct_coordinates[1:end-1,1].+duct_coordinates[2:end,1]), outputs.cp[1], label="FLOWFoil Duct with Center Body Effects",color=2) #hide
=#
```

Plotting the geometry and the output velocities and pressures show expected behavior when combining these two cases.

## Airfoil Cascade

For this example, we use data from chapter 2 of ["Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems" by  R. I. Lewis] (https://doi.org/10.1017/CBO9780511529542)

```@example cascade
using FLOWFoil

#this file contains the coordinates of the C4/70C50 airfoil as defined by lewis as well as the values of the pressure coefficient
data_path = normpath(
    joinpath(
        splitdir(pathof(FLOWFoil))[1], "..", "test", "data", "chapter2_lewis_validation.jl"
    ),
)
include(data_path)

#previously defined coordinates in chapter2_lewis_validation.jl
coordinates = [x y]

#setup Martensen method
method = Martensen(; solidity=1.0 / 0.900364, stagger=0.0)

#define flow angles (angles of attack) - in this case the angles of attack = the inflow angles hence why stagger is 0
flow_angles = [-35.0, 35.0]

#solve for outputs
outputs = analyze(coordinates, flow_angles; method=method)
```
