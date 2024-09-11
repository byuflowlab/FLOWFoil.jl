"""
    Mfoil{} <: Method
        
**Fields:**
- `viscous::Bool` : Flag whether to run an inviscid (false) or viscous (true) analysis.
"""
@kwdef struct Mfoil <: Method 
    viscous::Bool = false
end

const Xfoil = Mfoil
