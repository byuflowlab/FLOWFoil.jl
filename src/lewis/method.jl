"""
    Lewis{TB} <: Method

**Fields:**
- `body_of_revolution::Bool` : Flag(s) whether bodies are bodies of revolutions or not (ducts if not)
"""
@kwdef struct Lewis{TB} <: Method
    body_of_revolution::TB = [false]
end
