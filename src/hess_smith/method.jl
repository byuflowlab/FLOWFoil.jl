"""
    HessSmith <: Method

"""
@kwdef struct HessSmith{TF} <: Method 
    V_inf::TF=1.0
end