"""
    Lewis{TB} <: Method

# Fields:
- `body_of_revolution::Bool` : Flag(s) whether bodies are bodies of revolutions or not (`false` indicates an annular airfoil)

Note that if multiple bodies are used, the annular airfoil should come before the body of revolution.
"""
@kwdef struct Lewis{TB} <: Method
    body_of_revolution::TB = [false]
end
