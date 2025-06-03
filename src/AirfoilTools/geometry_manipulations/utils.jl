"""
    dot(A, B) = sum(a * b for (a, b) in zip(A, B))

A faster dot product.
"""
dot(A, B) = sum(a * b for (a, b) in zip(A, B))

"""
    norm(A) = sqrt(mapreduce(x -> x^2, +, A))

A faster 2-norm.
"""
norm(A) = sqrt(mapreduce(x -> x^2, +, A))

"""
    linear_transform(source_range, target_range, source_values)

Linear transfrom of values from range `(source_range[1], source_range[end])` to `(target_range[1], target_range[end])`

# Arguments
- `source_range::Vector{Float}` : range values come from (can also be a Tuple)
- `target_range::Vector{Float}` : range onto which we are transforming (can also be a Tuple)
- `source_values::Array{Float}` : array of source values to transform

# Returns
 - `target_values::Array{Float}` : array of transformed sourcevalues onto target range
"""
function linear_transform(source_range, target_range, source_values)
    return target_range[1] .+
           (target_range[end] - target_range[1]) .* (source_values .- source_range[1]) /
           (source_range[end] - source_range[1])
end
