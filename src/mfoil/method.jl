"""
    Mfoil <: Method

# Fields:
- `viscous::Bool` : Flag whether to run an inviscid (false) or viscous (true) analysis.
"""
@kwdef struct Mfoil{TB} <: Method
    viscous::TB = false
end

const Xfoil = Mfoil
