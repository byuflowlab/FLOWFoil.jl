"""
    Lewis <: Method

# Fields
- `body_of_revolution::AbstractVector{Bool}` : Flag(s) whether bodies are bodies of revolutions or not (`false` indicates an annular airfoil)

Note that if multiple bodies are used, the annular airfoil should come before the body of revolution.
"""
struct Lewis{TB} <: Method
    body_of_revolution::TB
end

function Lewis(; body_of_revolution=false)
    return Lewis(isscalar(body_of_revolution) ? [body_of_revolution] : body_of_revolution)
end
