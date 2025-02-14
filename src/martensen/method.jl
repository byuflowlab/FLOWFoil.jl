"""
    Martensen{} <: Method

# Fields:
- `pitch::Float` : < add description >
- `stagger::Float` : < add description >
"""
@kwdef struct Martensen{TP,TS} <: Method
    pitch::TP
    stagger::TS
end
