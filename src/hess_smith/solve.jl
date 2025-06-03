"""
    solve(::HessSmith, system_matrices)

Solves the linear system `A * γ = b` arising from the Hess-Smith panel method, where `A` is the
vortex influence matrix and `b` is the boundary condition vector.

# Arguments
- `::HessSmith`: The Hess-Smith solver type, used for dispatch.
- `system_matrices`: A named tuple with fields:
  - `A::Matrix`: The system matrix of vortex influence coefficients.
  - `b::Vector`: The boundary condition vector.

# Returns
- `γ::Vector`: The solution vector of vortex strengths (circulation values) at each panel.
"""
function solve(::HessSmith, system_matrices)
    return ImplicitAD.linear_solve(system_matrices.A, system_matrices.b)
end
