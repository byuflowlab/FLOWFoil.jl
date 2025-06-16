# Public API

## Analyze

```@docs; canonical=false
FLOWFoil.analyze
```

## Methods

FLOWFoil provides access to the following methods.

```@docs; canonical=false
FLOWFoil.Mfoil
FLOWFoil.Lewis
FLOWFoil.Martensen
FLOWFoil.HessSmith
FLOWFoil.LegacyXfoil
FLOWFoil.NeuralFoil
```

## Outputs

Based on the method used, FLOWFoil returns one of the following types:

```@docs; canonical=false
FLOWFoil.InviscidOutputs
FLOWFoil.LegacyXFOutputs
```

Or in the case of NeuralFoil, the native NeuralFoil output (`NeuralOutputs`) type will be returned
