# API Reference

# Public


--------------------------------------------


# Private

## System Types and Functions
```@docs
FLOWFoil.PlanarMesh
FLOWFoil.PlanarMeshSystem
FLOWFoil.AxiSymPanel
FLOWFoil.InviscidSystem
FLOWFoil.get_inviscid_system
FLOWFoil.size_system
FLOWFoil.assemble_boundary_conditions
FLOWFoil.assemble_ring_boundary_conditions
FLOWFoil.countkutta
```

## Influence Coefficient Functions
```@docs
FLOWFoil.assemble_vortex_matrix
FLOWFoil.assemble_vortex_coefficients
FLOWFoil.get_vortex_influence
FLOWFoil.get_source_influence
FLOWFoil.get_psibargamma
FLOWFoil.get_psitildegamma
FLOWFoil.get_psibarsigma
FLOWFoil.get_psitildesigma
FLOWFoil.assemble_ring_vortex_matrix
FLOWFoil.assemble_ring_vortex_coefficients
FLOWFoil.get_ring_vortex_influence
FLOWFoil.get_v_ring_vortex
FLOWFoil.get_u_ring_vortex
FLOWFoil.get_elliptics
FLOWFoil.get_offset
```

## Panel Geometries
```@docs
FLOWFoil.get_r
FLOWFoil.get_a
FLOWFoil.get_h
FLOWFoil.get_d
FLOWFoil.get_normal
FLOWFoil.get_theta
FLOWFoil.get_orientation
FLOWFoil.get_tangent
FLOWFoil.get_ring_geometry
FLOWFoil.get_distances
FLOWFoil.get_trailing_edge_info
```

## Solvers
```@docs
FLOWFoil.solve_inviscid
FLOWFoil.solve_inviscid_system
```

## Post Processing (2D)
```@docs
FLOWFoil.calculate_stream_grid
FLOWFoil.get_stream_grid_value
FLOWFoil.get_gamma_magnitudes
```

## Post Processing (Axisymmetric)
```@docs
FLOWFoil.get_relative_geometry_axisym
FLOWFoil.get_mesh_gammas
FLOWFoil.calculate_duct_thrust
```
