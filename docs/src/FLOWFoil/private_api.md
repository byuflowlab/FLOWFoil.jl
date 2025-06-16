# Private API

## Universal Dispatch

These functions are universally called in the convenience function, [analyze](@ref) except for the [NeuralFoil](@ref) and [LegacyXfoil](@ref) methods, which bypass these functions in favor of simply calling those packages' analysis functions.

```@docs
FLOWFoil.reformat_inputs
FLOWFoil.generate_panel_geometry
FLOWFoil.generate_system_geometry
FLOWFoil.generate_system_matrices
FLOWFoil.solve
FLOWFoil.post_process
```

## Universal Utility Functions

```@docs
FLOWFoil.get_d
FLOWFoil.get_r
FLOWFoil.get_panel_tangent
FLOWFoil.get_panel_normal
FLOWFoil.calculate_chord
FLOWFoil.linear_transform
FLOWFoil.smooth_distributions
FLOWFoil.dot
FLOWFoil.smooth_beta
FLOWFoil.laitone_compressibility_correction
```

## Additional Method Utilities

### Lewis
```@autodocs
Modules = [FLOWFoil]
Pages = ["lewis/geometry_utils.jl"]
```

### Mfoil
```@autodocs
Modules = [FLOWFoil]
Pages = ["mfoil/geometry_utils.jl"]
```
