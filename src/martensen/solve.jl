"""
    solve(::Martensen, system_matrices)

Solves the linear system for the Lewis method to find the strengths of singularities on panels.

# Arguments
- `::Martensen`: Marker type indicating the Martensen method.
- `system_matrices`: A NamedTuple containing:
  - `A`: The coefficient matrix (LHS).
  - `b`: The right-hand side vector or matrix.

# Returns
- `x`: The solution vector or matrix containing singularity strengths.
"""
function solve(::Martensen, system_matrices)
    return ImplicitAD.linear_solve(system_matrices.A, system_matrices.b)
end
