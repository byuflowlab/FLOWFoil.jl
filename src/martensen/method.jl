"""
    Martensen <: Method

# Fields
- `cascade::Bool=true` : flag to apply cascade treatment or not
- `solidity::Float=0.0` : Ratio between airfoil chord length and pitch. Airfoil pitch is simply the distance between chordlines in the cascade.
- `stagger::Float=0.0` : Angle (in degrees) from axis of rotation to airfoil chordline. Note that stagger is equivalent to the inflow angle minus the angle of attack.
- `transition_value::Float=Inf` : pitch to chord ratio at which we stop applying cascade effects (Lewis uses 30 in his implementation)
- `curvature_correction::Bool=false` : flag to apply curvature correction from Lewis
"""
@kwdef struct Martensen{TB,TF,TP,TS} <: Method
    cascade::TB = true
    solidity::TP = 0.0
    stagger::TS = 0.0
    transition_value::TF = Inf
    curvature_correction::TB = false
end
