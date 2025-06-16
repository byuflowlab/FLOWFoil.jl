# Public API

```@contents
Pages = ["public_api.md"]
Depth = 5
```

## Analyze

```@docs
FLOWFoil.analyze
```

## Methods

FLOWFoil provides access to the following methods.

### Hess Smith

```@docs
FLOWFoil.HessSmith
```

### LegacyXfoil

```@docs
FLOWFoil.LegacyXfoil
```

### Lewis

```@docs
FLOWFoil.Lewis
```

### Martensen

```@docs
FLOWFoil.Martensen
```

### Mfoil

```@docs
FLOWFoil.Mfoil
```

### NeuralFoil

```@docs
FLOWFoil.NeuralFoil
```

## Airfoil Tools

```@autodocs
module = [AirfoilTools]
order = [:function, :type]
```

## Outputs

Based on the method used, `FLOWFoil` returns one of the following types:

- `InvscidOutputs`: returned by [`HessSmith`](@ref), [`Lewis`](@ref), etc.
- `NeuralOutputs`: returned by [`NeuralFoil`](@ref)
- `LegacyXFOutputs`: returned by [`LegacyXfoil`](@ref)

```@docs
FLOWFoil.InvscidOutputs
FLOWFoil.NeuralOutputs
FLOWFoil.LegacyXFOutputs
```