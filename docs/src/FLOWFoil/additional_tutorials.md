# Tutorials

FLOWFoil includes various panel method implementations that are available based on the `method` keyword argument.
The analyze function is the way to run analyses with any method available in FLOWFoil, you just need to select the method you want and provide any additionally required inputs.

```@docs
FLOWFoil.analyze
```

## Xfoil Method

An Xfoil-like method, based on [Xfoil](https://websites.umich.edu/~kfid/codes.html), can be accessed using the `Xfoil` method type. It is also exported as `Mfoil` (an alias) since the implementation follows Drela's `mfoil` derivation. The method has two modes:

- **Inviscid** (default): a linear-vortex panel solve over the airfoil. Returns an `InviscidOutputs` object.
- **Viscous**: a coupled boundary-layer + inviscid solve. Set `viscous=true` and provide `Re`, `Mach`, and the other viscous knobs. Returns a NamedTuple (see [Viscous Mode](@ref) below).

### Inviscid Mode

```@example Xfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

method = Xfoil()

outputs = analyze(x, y, angles_of_attack; method=method)
```

The inviscid Xfoil outputs are of type `InviscidOutputs`. This is also the default method used in the [Quick Start](@ref).

### Viscous Mode

When `viscous=true` is set, `analyze` runs the coupled mfoil-style boundary-layer + inviscid panel solve described by [Drela](https://doi.org/10.1007/978-3-642-84010-4_1) and [Fidkowski](https://doi.org/10.2514/1.J061341). Each angle of attack is solved independently via an internal loop. A worked example with plots lives on the [Additional Examples](@ref) page.

```@example xfoil_viscous
using FLOWFoil

x, y = AirfoilTools.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)

method = Xfoil(viscous=true, Re=1e6, Mach=0.0, ncrit=9.0, etol=1e-10, maxiters=50)

outputs = analyze([x y], [5.0]; method=method)
nothing # hide
```

```@example xfoil_viscous
(cl=outputs.cl[1],
 cd=outputs.cd[1],
 cm=outputs.cm[1],
 transition_upper=outputs.transition_upper[1],
 transition_lower=outputs.transition_lower[1],
 converged=outputs.converged[1])
```

The viscous return is a flat NamedTuple. The top-level scalar fields per angle of attack are `cl`, `cd`, `cm`, `transition_upper`, `transition_lower`, and `converged`. The inviscid intermediates `gamma_ref` (node vortex strengths) and `A` (the linear-system influence matrix) are also exposed alongside. Per-station boundary-layer state for each α is in `outputs.per_alpha[i]` and includes `thetas_u/l/w`, `deltas_u/l/w`, `state3s_u/l/w` (laminar amplification factor `n` or turbulent shear-stress coefficient `cτ`), and the edge velocities `Veu`, `Vel`, `Vwake`.

!!! note "Viscous-mode constraints"
    - **Single body only.** A viscous run on a multi-body system raises an error.
    - **Blunt trailing edge required.** Use `AirfoilTools.naca4(...; blunt_trailing_edge=true)` or supply your own blunt-TE coordinates.

```@docs
FLOWFoil.Xfoil
```

## Lewis' Method for Axisymmetric Bodies

An axisymmetric method based on the one described by [Lewis](https://doi.org/10.1017/CBO9780511529542) can be accessed using the `Lewis` method type:

```@example lewis
using FLOWFoil

# Since this is an axisymmetric method, we'll use r instead of y
x, r = AirfoilTools.naca4()

# give the duct some diameter so it doesn't have negative radial dimensions (see warning below)
r .+= 1.0

# indicate that the body is not a body of revolution (i.e. a duct)
method = Lewis(; body_of_revolution=false)

# note: we need to input an an angle of attack, even though it is unused in this method.
outputs = analyze(x, r, 0.0; method=method)
```

The comments here mention multiple bodies, for more information, see the multi-body example: [Axisymmetric Mutli-element Systems](@ref) on the next page.

The outputs for the Lewis method are also of type InviscidOutputs.

```@docs
FLOWFoil.Lewis
```

!!! warning
    No part of the geometry for an axisymmetric body can reside below z=0, otherwise an error will be thrown.

## Martensen's Method for Axial Cascades

A periodic method for cascade analysis based on that developed by [Martensen](https://archive.org/details/nasa_techdoc_19710021012) can be accessed using the `Martensen` method type:

```@example martensen
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

# The cascade method requires solidity (closeness) of sections and stagger (inflow angle - angle of attack)
method = Martensen(solidity=1.2, stagger=15.0)

outputs = analyze(x, y, angles_of_attack; method=method)
```

The InviscidOutputs type is also used for the Martensen method.

```@docs
FLOWFoil.Martensen
```
!!! note
    If the `cascade` option is set to false, this method becomes a standard planar airfoil method, but uses constant vortices, so the Xfoil method is the superior method in that case.

## NeuralFoil Method

[NeuralFoil](https://github.com/peterdsharpe/NeuralFoil) is a multi-layer perceptron model of Xfoil.
We provide the Neuralfoil Method through the `NeuralFoil` method type and is accessed through the [NeuralFoil.jl](https://github.com/byuflowlab/NeuralFoil.jl) package:

```@example neuralfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

reynolds = 2e6
mach = 0.0

method = NeuralFoil(reynolds, mach; model_size="xlarge", n_crit=9.0, xtr_upper=1.0, xtr_lower=1.0)

outputs = analyze([x y], angles_of_attack; method=method)
```

```@docs
FLOWFoil.NeuralFoil
FLOWFoil.NeuralFoil(reynolds)
```

Note that the NeuralFoil method does not allow multi-body analysis like the other methods do as it is based specifically on Xfoil.  We also return a separate output type for the NeuralFoil method from the NeuralFoil.jl namespace.

## LegacyXfoil Method

We also have the LegacyXfoil method that is based on Xfoil and can be accessed with the `LegacyXfoil` method type:

```@example legacyxfoil
using FLOWFoil

x, y = AirfoilTools.naca4()

angles_of_attack = range(-5.0, 15.0; step=1)

reynolds = 2e6

method = LegacyXfoil(; reynolds=reynolds, npan=140)

outputs = analyze([x y], angles_of_attack; method=method)
```

Pass `reynolds=nothing` (or simply omit `reynolds`) to run the underlying Xfoil.jl in inviscid mode. When a single angle of attack is requested, `outputs.cp`/`outputs.x_cp` hold the surface pressure distribution and (in viscous mode) `outputs.bl` holds the boundary-layer dump `(s, x, y, ue, dstar, theta, cf)`; multi-α sweeps leave those fields as `nothing`.

```@docs
FLOWFoil.LegacyXfoil
FLOWFoil.LegacyXfoil()
```

Note that we return a separate output type for the LegacyXFoil method:

```@docs
FLOWFoil.LegacyXFOutputs
```

## Method comparison: NACA 2412 at Re=1e6, α=5°

The three Xfoil-flavored methods in FLOWFoil — `Xfoil(viscous=true)`, `LegacyXfoil` (wrapper around Xfoil.jl), and `NeuralFoil` (wrapper around NeuralFoil.jl) — can all return a pressure distribution and boundary-layer state for a single operating point. The example below runs all three on the same NACA 2412 / Re = 1e6 / α = 5° / M = 0 case and compares the surface cp and the momentum-thickness θ.

```@example xfoil_method_compare
include("../assets/plots_default.jl") #hide
using FLOWFoil
using Plots
using LaTeXStrings

x, y = AirfoilTools.naca4(2.0, 4.0, 12.0; N=201, blunt_trailing_edge=true)

out_vis = analyze([x y], [5.0]; method=Xfoil(viscous=true, Re=1e6, Mach=0.0,
                                              ncrit=9.0, etol=1e-10, maxiters=50))
out_lxf = analyze([x y], [5.0]; method=LegacyXfoil(; reynolds=1e6, npan=200))
out_nf  = analyze([x y], [5.0]; method=NeuralFoil(1e6; n_crit=9.0))
nothing # hide
```

Integrated coefficients:

```@example xfoil_method_compare
(
    Xfoil_viscous = (cl=out_vis.cl[1], cd=out_vis.cd[1], cm=out_vis.cm[1]),
    LegacyXfoil   = (cl=out_lxf.cl[1], cd=out_lxf.cd[1], cm=out_lxf.cm[1]),
    NeuralFoil    = (cl=out_nf.cl, cd=out_nf.cd, cm=out_nf.cm),
)
```

Surface pressure distribution. NeuralFoil does not output `cp` directly, so we derive it from its edge-velocity field via the incompressible relation `Cp = 1 − (Ue/V∞)²` (exact at M = 0). The NeuralFoil upper and lower samples are concatenated into a single curve traversing the airfoil from upper TE through LE to lower TE, matching the LegacyXfoil and FLOWFoil node ordering. Each method is drawn with a distinct linestyle so the (nearly-overlapping) curves remain visible:

```@example xfoil_method_compare
pa = out_vis.per_alpha[1]
cp_vis = pa.cp[1:length(x)]

cp_nf_upper = 1 .- vec(out_nf.upper_bl_ue_over_vinf).^2   # sampled LE → TE on upper
cp_nf_lower = 1 .- vec(out_nf.lower_bl_ue_over_vinf).^2   # sampled LE → TE on lower
x_nf_upper = range(0, 1; length=length(cp_nf_upper))
x_nf_lower = range(0, 1; length=length(cp_nf_lower))

# Stitch upper (reversed: TE → LE) + lower (LE → TE) into one continuous trace.
x_nf  = vcat(reverse(x_nf_upper),  x_nf_lower)
cp_nf = vcat(reverse(cp_nf_upper), cp_nf_lower)

plot(; xlabel=L"x/c", ylabel=L"c_p", yflip=true,
     title="NACA 2412, α=5°, Re=1e6", legend=:bottomright)
plot!(out_lxf.x_cp, out_lxf.cp; label="LegacyXfoil", linewidth=2)
plot!(x, cp_vis; label="FLOWFoil Xfoil (viscous)",
      linestyle=:dash, linewidth=2)
plot!(x_nf, cp_nf; label="NeuralFoil",
      linestyle=:dot, linewidth=2)
```

Boundary-layer momentum thickness θ. Each method's upper- and lower-surface traces both start at s = 0 (i.e. measured from the leading edge / stagnation point on each surface), so the airfoil region of the plot carries two sub-curves per method. The wake (for the two methods that produce one) is offset to start after the longer surface trace. Each method is drawn as a single Plots series with `NaN` breaks between sub-segments, giving one legend entry per method:

```@example xfoil_method_compare
# LegacyXfoil — split flat bl array at LE and at the wake junction (duplicated-TE row),
# then re-anchor each segment to s = 0.
bl = out_lxf.bl
i_wake = findfirst(iszero, diff(bl.s))
i_le   = argmin(view(bl.x, 1:i_wake))

lxf_s_upper = bl.s[i_le] .- reverse(bl.s[1:i_le])    # 0 (LE) → upper TE
lxf_t_upper = reverse(bl.theta[1:i_le])
lxf_s_lower = bl.s[i_le:i_wake] .- bl.s[i_le]        # 0 (LE) → lower TE
lxf_t_lower = bl.theta[i_le:i_wake]
wake_anchor_lxf = max(lxf_s_upper[end], lxf_s_lower[end])
lxf_s_wake = wake_anchor_lxf .+ (bl.s[(i_wake+1):end] .- bl.s[i_wake+1])
lxf_t_wake = bl.theta[(i_wake+1):end]

lxf_s = vcat(lxf_s_upper, NaN, lxf_s_lower, NaN, lxf_s_wake)
lxf_t = vcat(lxf_t_upper, NaN, lxf_t_lower, NaN, lxf_t_wake)

# FLOWFoil viscous — pa.su/pa.sl already start at 0 (stagnation); offset wake.
wake_anchor_vis = max(pa.su[end], pa.sl[end])
vis_s_wake = wake_anchor_vis .+ (pa.sw .- pa.sw[1])
vis_s = vcat(pa.su, NaN, pa.sl, NaN, vis_s_wake)
vis_t = vcat(pa.thetas_u, NaN, pa.thetas_l, NaN, pa.thetas_w)

# NeuralFoil — linearly map sample index onto each surface's arclength using LegacyXfoil's
# milestones. No wake.
nf_s_upper = range(0.0, lxf_s_upper[end]; length=length(cp_nf_upper))
nf_s_lower = range(0.0, lxf_s_lower[end]; length=length(cp_nf_lower))
nf_s = vcat(collect(nf_s_upper), NaN, collect(nf_s_lower))
nf_t = vcat(vec(out_nf.upper_theta), NaN, vec(out_nf.lower_theta))

plot(; xlabel=L"s", ylabel=L"\theta",
     title="NACA 2412, α=5°, Re=1e6", legend=:topleft)
plot!(lxf_s, lxf_t; label="LegacyXfoil", linewidth=2)
plot!(vis_s, vis_t; label="FLOWFoil Xfoil (viscous)",
      linestyle=:dash, linewidth=2)
plot!(nf_s, nf_t; label="NeuralFoil",
      linestyle=:dot, linewidth=2)
```

The viscous Xfoil and LegacyXfoil cp distributions agree closely (both solve the same coupled boundary-layer equations, just in different implementations); NeuralFoil's derived cp follows the same suction-peak shape with mild smoothing from its neural-network surrogate. With upper and lower θ both anchored at LE, the laminar-flat / transition-jump / post-transition-growth pattern lines up directly between methods on each surface.

