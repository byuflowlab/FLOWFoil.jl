"""
    Martensen{} <: Method

# Fields:
- `cascade::Bool` : flag to apply cascade treatment or not
- `pitch::Float` : < add description >
- `stagger::Float` : < add description >
- `transition_value::Float` : pitch_to_chord ratio at which we stop applying cascade effects
- `transition_hardness::Float` : FLOWMath.sigmoid_blend hardness for blend between planar and cascade influence coefficients
"""
@kwdef struct Martensen{TB,TF,TP,TS} <: Method
    cascade::TB = false
    pitch::TP = 0.0
    stagger::TS = 0.0
    transition_value::TF = 30.0
    transition_hardness::TF = 100.0
end
