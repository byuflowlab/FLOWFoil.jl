# Additional Examples

## Xfoil: Single inviscid airfoil comparision to analytic solution

This example uses the same Joukowsky airfoil presented in the original Xfoil paper.
We show here that our derivation and implementation of an Xfoil-like method also matches well to the analytical solution.

```@example Joukowsky
using FLOWFoil

center = [-0.1; 0.1]
radius = 1.0
alpha = 4.0

# - Joukowsky Geometry - #
x, y = FLOWFoil.AirfoilTools.joukowsky(center, radius, N=161)

# Plot Geometry
include("../assets/plots_default.jl") #hide
plot(x, y; aspectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

```@example Joukowsky
# - Analytic Solution - #
surface_velocity, surface_pressure_coefficient, cl = FLOWFoil.AirfoilTools.joukowsky_flow(
    center, radius, alpha; N=161
)

# - FLOWFoil Solution - #
outputs = analyze([x y], alpha; method=Xfoil())

# - Plot Outputs - #
include("../assets/plots_default.jl") #hide
plot(
    x[2:(end - 1)],
    surface_pressure_coefficient[2:(end - 1)];
    xlabel=L"x",
    ylabel=L"c_p",
    yflip=true,
    linestyle=:dash,
    linewidth=4, # hide
    label="Analytic Solution",
)
plot!(x[2:(end - 1)], outputs.cp[2:(end - 1)]; label="FLOWFoil")
```

---

## Axisymmetric Body of Revolution

For this example, we use experimental data of a body of revolution from chapter 4 of ["Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems" by  R. I. Lewis](https://doi.org/10.1017/CBO9780511529542).
Note that for the axisymmetric examples, we will use (z,r) rather than (x,y) to name our coordinates in the common cylindrical frame that ducts and bodies of revolution are defined.

```@example axisym
using FLOWFoil

data_path = normpath(
    joinpath(
        splitdir(pathof(FLOWFoil))[1], "..", "test", "data", "bodyofrevolutioncoords.jl"
    ),
)
# read in center_body_coordinates as well as normalized surface velocities.
include(data_path)

# - Plot Geometry - #
include("../assets/plots_default.jl") #hide
plot(
    center_body_coordinates[:, 1],
    center_body_coordinates[:, 2];
    aspectratio=1,
    ylim=(0,Inf),
    xlabel=L"z",
    ylabel=L"r",
    label=""
)
```

```@example axisym
# - FLOWFoil Solution - #
outputs = analyze(center_body_coordinates, [0.0]; method=Lewis(; body_of_revolution=[true]))

# - Plot Outputs - #
include("../assets/plots_default.jl") #hide
scatter(
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    label="Experimental Data",
    markerstrokecolor=1, # hide
    markercolor=1, # hide
    markersize=4, # hide
    xlabel=L"\frac{z}{c}",
    ylabel=L"\frac{V_s}{V_\infty}",
    legend=:bottomright, # hide
)
plot!(
    0.5 * (center_body_coordinates[1:(end - 1), 1] .+ center_body_coordinates[2:end, 1]),
    outputs.vs;
    label="FLOWFoil",
)
```

---

## Axisymmetric Annular Airfoil (Duct)

If we define an airfoil shape in an axisymmetric scheme, we model an annular airfoil, or in other words, a duct.  To do so, we follow a similar procedure to bodies of revolution with the exception that we set `body_of_revolution=false`.
The following case comes from the same book as the body of revolution above, but this time for a duct with a NACA 662-015 cross section.

```@example axisym
using FLOWFoil

duct_path = normpath(
    joinpath(splitdir(pathof(FLOWFoil))[1], "..", "test", "data", "naca_662-015.jl")
)
include(duct_path)

# - Plot Geometry - #
plot(
    duct_coordinates[:, 1],
    duct_coordinates[:, 2];
    aspectratio=1,
    xlabel=L"z",
    ylabel=L"r",
    label="",
)
```

```@example axisym
# - FLOWFoil Solution - #
outputs = analyze(duct_coordinates, [0.0]; method=Lewis(; body_of_revolution=[false]))

# - Plot Outputs - #
include("../assets/plots_default.jl") #hide
scatter(
    pressurexupper,
    pressureupper;
    markershape=:utriangle,
    label="Experimental Nacelle",
    color=1, # hide
    yflip=true,
    markerstrokecolor=1, # hide
    markercolor=1, # hide
    markersize=4, # hide
    xlabel=L"\frac{z}{c}",
    ylabel=L"c_p",
    legend=:bottomright, # hide
)
scatter!(
    pressurexlower,
    pressurelower;
    markershape=:dtriangle,
    label="Experimental Casing",
    color=1, # hide
    markerstrokecolor=1, # hide
    markercolor=1, # hide
    markersize=4, # hide
)
plot!(
    0.5 * (duct_coordinates[1:(end - 1), 1] .+ duct_coordinates[2:end, 1]),
    outputs.cp;
    label="FLOWFoil",
    color=2, # hide
)
```

As above, we plot experimental results along with our calculated values.


---

## Axisymmetric Mutli-element Systems

As an example of an multi-element axisymmetric system (such as that used for a ducted rotor), we will simply combine the two previous cases.
Note that we include the coordinates for the various bodies as a tuple of matrices, and in this case, we need to indicate in the `Lewis` method fields which of the bodies is a body of revolution (currently, the method only works if the annular airfoil comes first).
We simply put the previous two cases together to show that a multi-body cases works, and the results are reasonable.

```@example axisym
using FLOWFoil

# - Plot Geometry - #
plot(
    center_body_coordinates[:, 1],
    center_body_coordinates[:, 2];
    aspectratio=1,
    xlabel=L"z",
    ylabel=L"r",
    label="Center Body",
    legend=:right #hide
)
plot!(
    duct_coordinates[:, 1],
    duct_coordinates[:, 2];
    label="Duct",
)
```

```@example axisym
outputs = analyze(
    (duct_coordinates, center_body_coordinates), [0.0];
    method=Lewis(; body_of_revolution=[false, true]),
)

scatter(
    Vs_over_Vinf_x,
    Vs_over_Vinf_vs;
    label="Experimental Center Body",
    markerstrokecolor=1, # hide
    markercolor=1, # hide
    markersize=4, # hide
    xlabel=L"\frac{z}{c}",
    ylabel=L"\frac{V_s}{V_\infty}",
    legend=:bottomright, # hide
)
plot!(
    0.5 * (center_body_coordinates[1:(end - 1), 1] .+ center_body_coordinates[2:end, 1]),
    outputs.vs[2];
    label="FLOWFoil Center Body with Duct Effects",
)
```

```@example axisym
scatter(
    pressurexupper,
    pressureupper;
    markershape=:utriangle,
    label="Experimental Nacelle",
    color=1, #hide
    yflip=true,
    markerstrokecolor=1, #hide
    markercolor=1, #hide
    markersize=4, #hide
    xlabel=L"\frac{z}{c}",
    ylabel=L"c_p",
    legend=:bottomright, #hide
)
scatter!(
    pressurexlower,
    pressurelower;
    markershape=:dtriangle,
    label="Experimental Casing",
    color=1, #hide
    markerstrokecolor=1, #hide
    markercolor=1, #hide
    markersize=4, #hide
)
plot!(
    0.5 * (duct_coordinates[1:(end - 1), 1] .+ duct_coordinates[2:end, 1]),
    outputs.cp[1];
    label="FLOWFoil Duct with Center Body Effects",
    color=2, #hide
)
```

Plotting the geometry and the output velocities and pressures show expected behavior when combining these two cases.


---

## Linear Cascade

For this example, we uses data from chapter 2 of ["Vortex Element Methods for fluid Dynamic Analysis of Engineering Systems" by  R. I. Lewis] (https://doi.org/10.1017/CBO9780511529542) for a linear cascade section with inflow angle of +/- 35 degrees with zero twist.


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

# - Plot Geometry - #
include("../assets/plots_default.jl") #hide
plot(x, y; apsectratio=1, xlabel=L"x", ylabel=L"y", label="")
```

```@example cascade
# - Set up remaining inputs - #
# Define method inputs
method = Martensen(; solidity=1.0 / 0.900364, stagger=0.0)

# define angles of attack
# in this case the angles of attack = the inflow angles since stagger is 0
flow_angles = [-35.0, 35.0]

# - FLOWFoil Solution - #
outputs = analyze(coordinates, flow_angles; method=method)

# Panel midpoints (x has 51 nodes → 50 panels)
xmid = 0.5 .* (x[1:end-1] .+ x[2:end])

# - Plot Outputs - #
include("../assets/plots_default.jl") #hide
# Plot cp for alpha = -35 degrees

scatter(
    x_from_web_plot_digitizer_negative_35_degrees,
    cp_Lewis_negative_35_degrees;
    xlabel=L"\frac{x}{c}",
    ylabel=L"c_p",
    label="Lewis: -35°",
    markerstrokewidth=0, #hide
    yflip=true,
    markersize=4, #hide
    legend=:topright, #hide
)
plot!(xmid, outputs.cp[:, 1]; label="FLOWFoil")
```
```@example cascade
# Plot cp for alpha = +35 degrees

scatter(
    x_from_web_plot_digitizer_positive_35_degrees,
    cp_Lewis_positive_35_degrees;
    xlabel=L"\frac{x}{c}",
    ylabel=L"c_p",
    label="Lewis: +35°",
    markerstrokewidth=0, #hide
    yflip=true,
    markersize=4, #hide
    legend=:topright, #hide
)
plot!(xmid, outputs.cp[:, 2]; label="FLOWFoil")
```
