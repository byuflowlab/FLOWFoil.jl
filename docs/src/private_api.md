# Private API

```@contents
Pages = ["private_api.md"]
Depth = 5
```

## Universal Dispatch

```@docs
FLOWFoil.reformat_inputs
FLOWFoil.generate_panel_geometry
FLOWFoil.generate_system_geometry
FLOWFoil.generate_system_matrices
FLOWFoil.solve
FLOWFoil.post_process
```

## Universal Geometry Utilities

```@docs
FLOWFoil.get_d
FLOWFoil.get_r
FLOWFoil.get_panel_tangent
FLOWFoil.get_panel_normal
FLOWFoil.calculate_chord
```

## Universal Utilities

```@docs
FLOWFoil.linear_transform
FLOWFoil.smooth_distributions
FLOWFoil.dot
FLOWFoil.smooth_beta
FLOWFoil.laitone_compressibility_correction
```