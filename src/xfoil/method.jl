"""
    Xfoil <: Method

# Fields
- `viscous::Bool` : Flag whether to run an inviscid (false) or viscous (true) analysis.

NOTE: viscous method not yet implemented.
"""
@kwdef struct Xfoil{TB,TI,TF,TX,TZ} <: Method
    viscous::TB = false
    visualize::TB = false
    angle_of_attack::TF = 0.0
    x_range::AbstractVector{TX} = [-2, 2]
    y_range::AbstractVector{TZ} = [-2, 2]
    Nx::TI = 100
    Ny::TI = 100
end

const Mfoil = Xfoil
