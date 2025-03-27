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


