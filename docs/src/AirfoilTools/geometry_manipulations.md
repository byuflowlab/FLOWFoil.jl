# Airfoil Geometry Manipulation Tools

Here we include the variety of methods for manipulating airfoil geometries in useful ways implemented in AirfoilTools.

## Deconstruction

```@docs
AirfoilTools.split_upper_lower
```

## Transformation

```@docs
AirfoilTools.flip!
AirfoilTools.zero_z_te!
AirfoilTools.rotate_coordinates!
AirfoilTools.normalize_coordinates!
AirfoilTools.position_coordinates!
```

## Re-definition

```@docs
AirfoilTools.whole_cosine_spacing
AirfoilTools.split_cosine_spacing
AirfoilTools.repanel_airfoil
AirfoilTools.refine_trailing_edge
```

## Contributing

We welcome the addition of more convenience functions for airfoil geometry manipulation.
