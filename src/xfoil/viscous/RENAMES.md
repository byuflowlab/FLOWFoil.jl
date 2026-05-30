# Viscous module renames (vs XfoilTyler)

When the XfoilTyler viscous code was ported into FLOWFoil, several function
and parameter names were cleaned up. This document lists what changed.

The conventions applied:
- Drop the `get_` prefix on accessor-style functions (Julia idiom).
- Replace XfoilTyler-legacy or abbreviated names with descriptive ones
  where the legacy name was unclear.
- Keep `get_` on closures whose returned value's natural local-variable
  name matches the function name — dropping `get_` there would cause
  `Hk = Hk(...)`-style shadowing that Julia rejects.
- Drop "Java-y" `make_` / `build_` prefixes on functions that just return
  a value (NamedTuple constructors).

All XfoilTyler-side mathematical symbols (`δ*`, `θ`, `ue`, `Hk`, `H*`,
`Reθ`, `Cτ0`, `Eτ0`, `KTb`, `KTl`, etc.) are kept because anyone reading
this code will be reading Drela's documents alongside.

## Function renames

### Geometry / wake (viscous_geometry.jl, viscous_wake.jl)

| XfoilTyler | FLOWFoil |
|---|---|
| `parseairfoil` | `find_stagnation_split` |
| `split_s` | `resplit_at_stagnation` |
| `get_s_length_vectors` | `arclength_vectors` |
| `geometric_spaced_points` | `geometric_spacing` |
| `get_TE_panel` | `build_te_geometry` |
| `get_panels` | `build_viscous_airfoil_panels` |

### Closures (viscous_closures.jl)

| XfoilTyler | FLOWFoil |
|---|---|
| `get_H` | `H_shape` |
| `get_δ` | `bl_thickness` |
| `get_Hstarstar` | `Hstarstar` |
| `get_Hstar_laminar` | `Hstar_laminar` |
| `get_Hstar_turbulent` | `Hstar_turbulent` |
| `get_cf_laminar` | `cf_laminar` |
| `get_cf_turbulent` | `cf_turbulent` |
| `get_cf_wake` | `cf_wake` |
| `get_cdi_laminar` | `cdi_laminar` |
| `get_cdi_turbulent` | `cdi_turbulent` |
| `get_cdi_lamwake` | `cdi_lamwake` |
| `get_cdi_wall` | `cdi_wall` |
| `get_cdi_outer` | `cdi_outer` |
| `get_cdi_stress` | `cdi_stress` |
| `get_ct_transition` | `ct_at_transition` |
| `get_dndx` | `amplification_rate` |
| `get_ηup` | `upwind_factor` |
| `get_uq` | `uq` |
| `apply_upwinding` | `upwind_avg` |

### Residuals (viscous_residuals.jl)

| XfoilTyler | FLOWFoil |
|---|---|
| `ResidualParamsStation` struct + functor | `bl_station(...)` returning a NamedTuple + `station_residual(states, station)` |
| `residuals_together_local` | `bl_station_residuals` |
| `residuals_transition` | `residuals_transition` (kept) |
| `residual_momentum` | `residual_momentum` (kept) |
| `residual_shape` | `residual_shape` (kept) |
| `residual_amplification` | `residual_amplification` (kept) |
| `residual_shearlag` | `residual_shearlag` (kept) |
| `get_residual_velocity` | `velocity_residual` |
| `solve_for_xt_implicit` / `xt_residual_implicit` | `_solve_for_xt_implicit` / `_xt_residual_implicit` (now private) |
| `xt_residual` | `_xt_residual_scalar` (private) |
| `local_newton_xt` | `_newton_for_xt` (private) |

### Solvers & initialization (viscous_initialization.jl, viscous_solvers.jl)

| XfoilTyler | FLOWFoil |
|---|---|
| `make_station` / `ResidualParamsStation(...)` | `bl_station(...)` |
| `local_newton` | `local_newton` (kept; now takes the station NamedTuple) |
| `station_laminar` | `station_laminar` (kept) |
| `station_turbulent` | `station_turbulent` (kept) |
| `station_stagnation` | `station_stagnation` (kept) |
| `station_transition` | `station_transition` (kept) |
| `initialize_boundary_layer` | `initialize_boundary_layer` (kept) |
| `initialize_wake_first_node` | `initialize_wake_first_node` (kept) |
| `interpolate_stagnation_state` | `interpolate_stagnation_state` (kept) |
| `GlobalResiduals` struct + functor | `build_global_residuals(...)` returning a NamedTuple + `global_residual(U, gr)` |
| `global_newton` | `global_newton` (kept) |
| `stagnation_point_move!` | `stagnation_point_move!` (kept) |
| `resolve_transition_location!` | `resolve_transition_location!` (kept) |
| `march_amplification!` | `march_amplification!` (kept) |
| `ResolveTransitionParams` struct + functor + `local_newton_amp` | inlined as `_amplification_residual` + `local_newton_amp` |

### Coupling (viscous_coupling.jl)

| XfoilTyler | FLOWFoil |
|---|---|
| `run_viscous` (single body) | `analyze_viscous(method, …)` + `_solve_one_alpha` (per-α loop) |
| `run_coupled_system_global` (multi entry-point) | dispatch via `analyze(coords, αs; method=Xfoil(viscous=true, …))` |
| `get_cl`, `get_cd`, `get_cdp` | `_viscous_cl`, `_viscous_cd`, `_viscous_cm` (private helpers) |

### mfoil influence integrals (viscous_influence.jl)

| XfoilTyler | FLOWFoil |
|---|---|
| `get_beta_mfoil` | `mfoil_beta` |
| `get_relative_geometry_panel_point` | `mfoil_relative_geometry` |
| `get_I1` (linear vortex Ψ) | inlined inside `mfoil_linear_vortex_psi` |
| `get_I3` (constant source Ψ) | inlined inside `mfoil_constant_source_psi` |
| `get_constant_source_streamfunction_influence` | `mfoil_constant_source_psi` |
| `get_constant_source_velocity_influence` | `mfoil_constant_source_velocity` |
| `get_linear_vortex_streamfunction_influence` | `mfoil_linear_vortex_psi` |
| `get_linear_vortex_velocity_influence` | `mfoil_linear_vortex_velocity` |
| `get_linear_source_streamfunction_influence` | `mfoil_linear_source_psi` |
| `get_linear_source_velocity_influence` | `mfoil_linear_source_velocity` |

## Functions that *kept* the `get_` prefix

Dropping `get_` from these would force `name = name(...)` patterns (e.g.,
`Hk = Hk(H, M2)`) all over the codebase. Julia treats the LHS as a fresh
local binding, which then shadows the global function on the RHS and
errors with `UndefVarError`. Renaming the local variables instead would
mean spreading suffixes like `Hk_val` through several hundred lines of
math code, which is less readable than keeping a `get_` on the function.

The closures below keep `get_`:

- `get_uk` — KT-corrected edge velocity. Result stored as `uk`.
- `get_Mach2` — local Mach² (result stored as `M2`, but other callers use `Mach2 = …`).
- `get_density` — local density. Result stored as `rho` or `density`.
- `get_Reθ` — momentum-thickness Reynolds. Result stored as `Reθ`.
- `get_Hk`, `get_Hkc` — kinematic shape parameter and its corrected form.
- `get_Us` — normalized wall slip velocity.
- `get_cτeq` — equilibrium √cτ for the shear-lag equation.
- Amplification rate building blocks: `get_Ĥ`, `get_L0`, `get_sn`,
  `get_rn`, `get_fn`, `get_gn`, `get_ϵn`. These are all nested inside
  `amplification_rate` where each is a `name = get_name(prev_name)` step.

## Local variable renames

In `cdi_turbulent(Hk, Reθ, cf, Us, Hstar, cτ; airfoil=true)` the per-component
locals were renamed to avoid shadowing the helper functions of the same
names:

| Before | After |
|---|---|
| `cdi_wall = cdi_wall(…)`     | `wall = cdi_wall(…)` |
| `cdi_outer = cdi_outer(…)`   | `outer = cdi_outer(…)` |
| `cdi_stress = cdi_stress(…)` | `stress = cdi_stress(…)` |
| `cdi_laminar = cdi_laminar(…)` | `laminar = cdi_laminar(…)` |

In `residual_shearlag`:

| Before | After |
|---|---|
| `uq = uq(Hk, cf, Reθ, δstar_avg; airfoil)` | `uq_val = uq(Hk, cf, Reθ, δstar_avg; airfoil)` |

## State / parameter naming

Some structural names changed in the port:

- XfoilTyler's mutable struct `OperatingParameters` was replaced by an
  immutable NamedTuple `thermo` (built by `init_viscous_thermo`) plus a
  small NamedTuple `state_refs` holding `Ref`s for the per-alpha mutable
  scalars (`sft_l`, `sft_u`, `transition_lower`, `transition_upper`,
  `wake_gap`).
- XfoilTyler's `Panel` / `WakePanel` / `TEgeometry` structs were
  replaced by NamedTuples returned from `build_viscous_airfoil_panels`,
  `build_wake`, and `build_te_geometry` — same field names so the math
  ports cleanly.
- The third BL state at each node is still called `state3` (either `n`
  in laminar flow or `cτ` in turbulent flow), matching XfoilTyler.
