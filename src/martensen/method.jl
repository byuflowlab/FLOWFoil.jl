"""
    Martensen{} <: Method

# Fields:

"""
@kwdef struct Martensen{TP,TS} <: Method
    pitch::TP
    stagger::TS
end
