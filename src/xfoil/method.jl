"""
    Xfoil <: Method

# Fields
- `viscous::Bool` : Flag whether to run an inviscid (false) or viscous (true) analysis.

NOTE: viscous method not yet implemented.
"""
@kwdef struct Xfoil{TB} <: Method
    viscous::TB = false
end

const Mfoil = Xfoil