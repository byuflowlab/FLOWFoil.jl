"""
    solve(::Lewis, system_matrices)

Solves the linear system for the Lewis method to find the strengths of singularities on panels.

# Arguments
- `::Lewis`: Marker type indicating the Lewis method.
- `system_matrices`: A NamedTuple containing:
  - `A`: The coefficient matrix (LHS).
  - `b`: The right-hand side vector or matrix.

# Returns
- `x`: The solution vector or matrix containing singularity strengths.
"""
function solve(::Lewis, system_matrices)
    #=
    # # Solve System
    if ndims(system_matrices.b) > 1
        # initialize output
        x = similar(system_matrices.b)
        
        # Pre-compute factorization
        Afact = LinearAlgebra.factorize(ImplicitAD.fd_value(system_matrices.A))
        # loop through dimensions, reusing the factorization of A
        for i in 1:size(system_matrices.b, 2)
            x[:, i] = ImplicitAD.implicit_linear(
            system_matrices.A, system_matrices.b[:, i]; Af=Afact
            )
        end
    else
        # Just solve like normal
        x = ImplicitAD.implicit_linear(system_matrices.A, system_matrices.b)
    end

    return x
    =#
    return ImplicitAD.linear_solve(system_matrices.A, system_matrices.b)
end
