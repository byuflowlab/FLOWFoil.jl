"""
    HessSmith <: Method

# Fields
- `V_inf::TF=1.0` : magnitude of the free stream velocity
"""
@kwdef struct HessSmith{TF} <: Method 
    V_inf::TF=1.0
end