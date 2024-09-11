#=
General Post Processing Types and Functions

Authors: Judd Mehr,
=#

abstract type Polar end

"""
    post_process(::ProblemType, problem, panels, mesh, solution; debug=false)

Post-process solution and produce a Polar object.

**Arguments:**
- `method::ProblemType` : Problem type for dispatch
- `problem::Problem` : Problem object
- `panels::Vector{Panel}` : vector of Panel objects
- `mesh::Mesh` : Mesh object
- `solution::Solution` : Solution object

**Keyword Arguments:**
- `npanels::Int` : number of panels to use on top and bottom surface for surface distribution smoothing
- `debug::Bool` : Flag for output format (see below)

**Returns:**
- If debug == false: `polar::Polar` : a Polar object
- If debug == true: all the fields of a Polar object (in order), but as a tuple rather than a struct, such that the raw, unsmoothed, velocity and pressure distributions can be output.
"""
function post_process(
    ::ProblemType, problem, panels, mesh, solution; npanels=80, debug=false
) end

######################################################################
#                                                                    #
#                       PLANAR POST PROCESSING                       #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#


#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

######################################################################
#                                                                    #
#                    AXISYMMETRIC POST PROCESSING                    #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    AxisymmetricPolar{TF,TA}

**Fields:**
- `surface_velocity::Matrix{Float}` : surface velocity on each panel
- `surface_pressure::Matrix{Float}` : surface pressure coefficient on each panel
- `xsmooth::Matrix{Float}` : x-values associated with smoothed surface distributions
"""
struct AxisymmetricPolar{TF} <: Polar
    surface_velocity::Matrix{TF}
    surface_pressure::Matrix{TF}
    xsmooth::Matrix{TF}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

function post_process(
    ap::AxisymmetricProblem, problem, panels::TP, mesh, solution; npanels=80, debug=false
) where {TP<:Panel}
    return post_process(
        ap::AxisymmetricProblem, problem, [panels], mesh, solution; npanels=80, debug=false
    )
end
function post_process(
    ap::AxisymmetricProblem, problem, panels, mesh, solution; npanels=80, debug=false
)

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    nbodies = mesh.nbodies

    # - Initialize Outputs - #
    TF = eltype(mesh.m)
    v_surf = zeros(TF, nbodies, 2 * npanels - 1)
    p_surf = zeros(TF, nbodies, 2 * npanels - 1)
    xsmooth = zeros(TF, nbodies, 2 * npanels - 1)

    for m in 1:nbodies
        # - Extract surface velocity - #
        #= Note that we here assume that we are using the subtractive method for the kutta conditions, requiring us to recover the solution value for the last panel (trailing edge upper side).  We also assume here that the indexing starts at the lower side trailing edge and proceeds clockwise back to the upper side trailing edge.  Otherwise, not only will the solver not have worked, there will also be an indexing error here
        =#
        vti = [solution.x[idx[m][1:(end - 1)]]; -solution.x[idx[m][1]]]

        # - Calculate surface pressure - #
        cpi = 1.0 .- (vti) .^ 2

        ### --- Smooth Distributions --- ###
        #smooth_distributions functions are found in utils.jl
        v_surf[m, :], xsmooth[m, :] = smooth_distributions(
            Constant(),
            panels[m].panel_center,
            vti,
            npanels;
            body_of_revolution=ap.body_of_revolution[m],
        )
        p_surf[m, :], _ = smooth_distributions(
            Constant(),
            panels[m].panel_center,
            cpi,
            npanels;
            body_of_revolution=ap.body_of_revolution[m],
        )
    end

    return AxisymmetricPolar(v_surf, p_surf, xsmooth)
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
function calculate_duct_thrust(polar::AxisymmetricPolar, panels, mesh; Vinf=1.0, rho=1.225)

    # - Rename for Convenience - #
    idx = mesh.node_indices

    # Calculate dynamic pressure
    q = 0.5 * rho * Vinf^2

    # Initialize output
    TF = eltype(mesh.m)
    fx = zeros(TF, mesh.nbodies)

    # Loop through Bodies
    for m in 1:(mesh.nbodies)
        if !problem.body_of_revolution[m]

            #dimensionalize pressure
            P = polar.surface_pressure[idx[m]] .* q

            # add panel pressure in x-direction to total sectional force
            fx += sum(
                P .* panels[m].panel_length[:] .* panels[m].panel_normal[:, 1] .*
                panels[m].center_point[:, 2],
            )
        end
    end

    #return total duct thrust for whole annulus: -fx*2pi
    return fx * 2.0 * pi
end
######################################################################
#                                                                    #
#                     PERIODIC POST PROCESSING                       #
#                                                                    #
######################################################################

#---------------------------------#
#              TYPES              #
#---------------------------------#

"""
    PeriodicPolar{TF,TA}

**Fields:**
- `surface_velocity::Matrix{Float}` : surface velocity on each panel
- `surface_pressure::Matrix{Float}` : surface pressure coefficient on each panel
- `xsmooth::Matrix{Float}` : x-values associated with smoothed surface distributions
"""
struct PeriodicPolar{TF} <: Polar
    surface_velocity::Array{TF,3}
    surface_pressure::Array{TF,3}
    xsmooth::Array{TF,3}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

function post_process(
    ap::PeriodicProblem, problem, panels::TP, mesh, solution; npanels=80, debug=false
) where {TP<:Panel}
    return post_process(
        ap::PeriodicProblem, problem, [panels], mesh, solution; npanels=80, debug=false
    )
end

function post_process(
    ::PeriodicProblem, problem, panels, mesh, solution; npanels=80, debug=false
)

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    nbodies = mesh.nbodies
    flow_angle = problem.flow_angle
    naoa = length(flow_angle)

    # - Initialize Outputs - #
    TF = eltype(mesh.x)
    # Surface Velocities
    v_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    # Surface Pressures
    p_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    xsmooth = zeros(TF, nbodies, 2 * npanels - 1, naoa)

    # vortex strengths
    gamma0 = solution.x[:, 1]
    gamma90 = solution.x[:, 2]

    # Flow angles
    flow_angle = problem.flow_angle

    for m in 1:nbodies
        for a in 1:naoa
            # - Extract surface velocity - #
            vti = solution.x[1:idx[end][end]]
            # vti = [
            #     gamma0[i] * cosd(flow_angle[a]) + gamma90[i] * sind(flow_angle[a]) for
            #     i in idx[m]
            # ]

            # - Calculate surface pressure - #
            cpi = 1.0 .- (vti) .^ 2

            ### --- Smooth Distributions --- ###
            #smooth_distributions functions are found in utils.jl
            #smooth_distributions functions are found in utils.jl
            v_surf[m, :, a], xsmooth[m, :, a] = smooth_distributions(
                Constant(), panels[m].panel_center, vti, npanels
            )
            p_surf[m, :, a], _ = smooth_distributions(
                Constant(), panels[m].panel_center, cpi, npanels
            )
        end
    end

    return PeriodicPolar(v_surf, p_surf, xsmooth)
end
