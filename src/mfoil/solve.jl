function analyze(method::Mfoil, system)

    if method.viscous
        error("Viscous method not yet implemented.")
    else
        return analyze_inviscid(method, system)
    end

end

"""
    analyze_inviscid(inviscid_system, mesh)

analyze the inviscid_system for the vortex and streamfunction strengths.

Outputs the InviscidSolution object which contains the inviscid_system.

**Arguments:**
- `inviscid_system::InviscidSystem` : inviscid_system to analyze.
- `mesh::Mesh` : Mesh defining geometry (to put into solution object)

**Returns:**
 - `x::inviscid solution, consisting of node vortex strengths`

"""
function analyze_inviscid(method::Mfoil, inviscid_system)

    # analyze System
    if ndims(inviscid_system.b) > 1

        # initialize output
        x = similar(inviscid_system.b)

        # Pre-compute factorization
        Afact = LinearAlgebra.factorize(ImplicitAD.fd_value(inviscid_system.A))

        # loop through dimensions, reusing the factorization of A
        for i in 1:size(inviscid_system.b)[2]
            x[:, i] = ImplicitAD.implicit_linear(
                inviscid_system.A, inviscid_system.b[:, i]; Af=Afact
            )
        end

    else
        # Just solve like normal
        x = ImplicitAD.implicit_linear(inviscid_system.A, inviscid_system.b)
    end

    return x
end
