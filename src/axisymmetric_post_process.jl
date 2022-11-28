#=
Axisymmetric Panel Method Post Processing

Authors: Judd Mehr,

Date Started: 04 November 2022

Change Log:
=#

"""
    AxiSymPolar{TF,TA}

**Fields:**
- `thrust::Float` : Thrust (or drag) of body
- `surface_velocity::Array{Float}` : surface velocity on each panel
- `surface_pressure::Array{Float}` : surface pressure coefficient on each panel
"""
struct AxiSymPolar{TF,TA} <: Polar
    thrust::TF
    surface_velocity::TA
    surface_pressure::TA
end

"""
    get_axisymmetric_polar(inviscid_solution, Vinf; rho=1.225)

Assemble post processing values for axisymmetric case.

**Arguments:**
- `inviscid_solution::FLOWFoil.InviscidSolution` : Inviscid solution object from FLOWFoil.
- `Vinf::Float` : freestream velocity

**Keyword Arguments:**
- `rho::Float` : air density

**Returns:**
- `axisym_post::AxiSymPost` : Post processed object for Axisymmetric cases.
"""
function get_axisymmetric_polar(inviscid_solution; Vinf=1.0, rho=1.225)
    surface_velocity = inviscid_solution.panelgammas
    surface_pressure = 1.0 .- (surface_velocity) .^ 2
    if inviscid_solution.meshes[1].bodyofrevolution
        thrust = calculate_duct_thrust(inviscid_solution; Vinf=Vinf, rho=1.225)
    else
        TF = eltype(surface_velocity)
        thrust = TF(0.0)
    end

    return AxiSymPolar(thrust, surface_velocity, surface_pressure)
end

"""
    calculate_duct_thrust(inviscid_solution; Vinf=1.0, rho=1.225)

Calculate the thrust of the duct.
TODO: need to test!

**Arguments:**
- `inviscid_solution::FLOWFoil.InviscidSolution` : Inviscid solution object from FLOWFoil.
- `Vinf::Float` : freestream velocity

**Keyword Arguments:**
- `rho::Float` : air density

**Returns:**
- `thrust::Float` : duct thrust (negative indicates drag)
"""
function calculate_duct_thrust(inviscid_solution; Vinf=1.0, rho=1.225)

    #unpack for convenience
    gammas = inviscid_solution.panelgammas
    duct_mesh = inviscid_solution.meshes[1]

    #calculate dynamic pressure
    q = 0.5 * rho * Vinf^2

    #initialize output
    fx = 0.0

    #get gammas specific to duct mesh (assumed to be first mesh)
    #note that FLOWFoil solves things in terms of 1/Vinf, so these surface velocities are actually Vti/Vinf already.
    duct_surface_velocities = get_mesh_gammas(gammas, inviscid_solution.meshes, 1)

    # loop through panels for this mesh
    for j in 1:length(duct_mesh.panels)

        #get current panel
        panel = duct_mesh.panels[j]

        #calculate pressure on panel
        cp_panel = 1.0 - (duct_surface_velocities[j])^2

        #dimensionalize
        P = cp_panel * q

        #add panel pressure in x-direction to total sectional force
        fx += P * panel.length * panel.normal[1] * panel.controlpoint[2]
    end

    #return total duct thrust for whole annulus: -fx*2pi
    return fx * 2.0 * pi
end

"""
    get_mesh_gammas(gammas, meshes, meshidx)

Get the gamma values only for the mesh at index meshidx in meshes.

**Arguments:**
- `gammas::FLOWFoil.InviscidSolution.panelgammas` : vortex strengths at each panel in the system.
- `meshes::Array{FLOWFoil.AxiSymMesh}` : Array of meshes in system
- `meshidx::Int` : index of which mesh in the meshes array for which to obtain the associated gammas.

**Returns:**
- `mesh_gammas::Array{Float}` : panel gamma values for input mesh
"""
function get_mesh_gammas(gammas, meshes, meshidx)

    #initialize offset
    offset = 0

    #if we're interested in values on mesh greater than 1, add to offset
    if meshidx > 1
        for i in 1:(meshidx - 1)
            offset += length(meshes[i].panels)
        end
    end

    #grab the gammas for just the body we want.
    mesh_gammas = gammas[(1 + offset):(offset + length(meshes[meshidx].panels))]

    return mesh_gammas
end
