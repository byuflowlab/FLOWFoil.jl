# API Reference

Additional Types and Functions not already shown in other parts of the documentation.

## Public

### Geometry and Meshing

#### General
```@docs
FLOWFoil.position_coordinates
FLOWFoil.position_coordinates!
```

#### 2D (Planar) systems
```@docs
FLOWFoil.PlanarMesh
FLOWFoil.PlanarMeshSystem
```

#### Axisymmetric Systems
```@docs
FLOWFoil.AxiSymPanel
```

### Problems and Solutions
```@docs
FLOWFoil.Problem
FLOWFoil.InviscidSolution
```

### Post Processing
```@docs
FLOWFoil.Polar
FLOWFoil.calculate_stream_grid
```

### Convenience Functions

#### Conformal Mapped Airfoils
```@docs
FLOWFoil.joukowsky
FLOWFoil.joukowskyflow
FLOWFoil.karman_trefftz
```

#### NACA 4-Series
```@docs
FLOWFoil.naca4
```

--------------------------------------------


## Private

### Inviscid System Construction

#### Setup
```@docs
FLOWFoil.size_system
FLOWFoil.get_offset
FLOWFoil.get_trailing_edge_info
FLOWFoil.countkutta
```

#### Panel Geometry
```@docs
FLOWFoil.get_r
FLOWFoil.get_d
FLOWFoil.get_theta
FLOWFoil.get_h
FLOWFoil.get_a
FLOWFoil.get_tangent
FLOWFoil.get_normal
FLOWFoil.get_distances
FLOWFoil.get_orientation
FLOWFoil.get_ring_geometry
```

#### Singularity Influences
```@docs
FLOWFoil.get_psibargamma
FLOWFoil.get_psitildegamma
FLOWFoil.get_vortex_influence
FLOWFoil.get_psibarsigma
FLOWFoil.get_psitildesigma
FLOWFoil.get_source_influence
FLOWFoil.get_ring_vortex_influence
```

#### Matrix Assembly

##### General
```@docs
FLOWFoil.InviscidSystem
FLOWFoil.get_inviscid_system
```

##### 2D (Planar)
```@docs
FLOWFoil.assemble_vortex_coefficients
FLOWFoil.assemble_vortex_matrix
FLOWFoil.assemble_boundary_conditions
```

##### Axisymmetric
```@docs
FLOWFoil.assemble_ring_vortex_coefficients
FLOWFoil.get_u_ring
FLOWFoil.get_v_ring
FLOWFoil.get_elliptics
FLOWFoil.assemble_ring_vortex_matrix
FLOWFoil.assemble_ring_boundary_conditions
```

### Inviscid System Solving
```@docs
FLOWFoil.solve_inviscid
FLOWFoil.solve_inviscid_system
```

### Post Processing
```@docs
FLOWFoil.get_vortex_magnitudes
FLOWFoil.get_stream_grid_value
FLOWFoil.get_gamma_magnitudes
```

### Convenience Functions
```@docs
FLOWFoil.cosinespacing
FLOWFoil.split_upper_lower
FLOWFoil.normalize_airfoil!
FLOWFoil.thickness
FLOWFoil.camber
```
