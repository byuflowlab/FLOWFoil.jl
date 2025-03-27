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

"""
    PlanarPolar{TF}

Also used for Periodic (cascade) post processing.

**Fields:**
- `lift::Matrix{Float}` : Lift Coefficient.
- `drag::Matrix{Float}` : Total Drag Coefficient.
- `pdrag::Matrix{Float}` : Pressure Drag Coefficient.
- `idrag::Matrix{Float}` : Induced Drag Coefficient.
- `moment::Matrix{Float}` : Moment Coefficient.
- `surfacevelocity::Array{Float}` : smoothed surface velocity distribution
- `surfacepressure::Array{Float}` : smoothed surface pressure distribution
- `xsmooth::Array{Float}` : x-values associated with smoothed surface distributions
"""
struct PlanarPolar{TF} <: Polar
    lift::Matrix{TF}
    drag::Matrix{TF}
    pdrag::Matrix{TF}
    idrag::Matrix{TF}
    moment::Matrix{TF}
    surface_velocity::Array{TF,3}
    surface_pressure::Array{TF,3}
    xsmooth::Array{TF,3}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

#TODO: once again, this is xfoil specific. Probably want to move this elsewhere eventually
function post_process(
    ::PlanarProblem, problem, panels, mesh, solution; npanels=80, debug=false
)

    ##### ----- Set Up ----- #####

    ### --- Rename for Convenience --- ###

    # number of airfoils
    nbodies = mesh.nbodies

    # number of angles of attack
    flow_angle = problem.flow_angle
    naoa = length(flow_angle)

    # node indices
    pidx = mesh.panel_indices
    nidx = mesh.node_indices

    # vortex strengths
    gamma0 = solution.x[:, 1]
    gamma90 = solution.x[:, 2]

    # chord length
    chord = mesh.chord

    ### --- Initialize Outputs --- ###

    # output floating point type
    TF = typeof(chord)

    # - Coefficients - #

    # Lift coefficient
    cl = zeros(TF, nbodies, naoa)

    # Total drag coefficient
    cd = zeros(TF, nbodies, naoa)
    # Profile drag coefficient
    cdp = zeros(TF, nbodies, naoa)
    # Inviscid drag coefficient
    cdi = zeros(TF, nbodies, naoa)

    # Moment coefficient
    cm = zeros(TF, nbodies, naoa)

    # - Surface Distributions - #

    #= NOTE:
        If debug is true, then we can allow for the raw outputs for surface distributions which may or may not be the same length (thus requiring a tuple).
        If debug is false, then we assume that we'll be smoothing the outputs with Akima splines.
    =#
    if debug
        # Surface Velocities
        v_surf = NTuple(nbodies, Matrix{TF}(undef, (2 * npanels - 1, naoa)))
        # Surface Pressures
        p_surf = NTuple(nbodies, Matrix{TF}(undef, (2 * npanels - 1, naoa)))
    else
        # Surface Velocities
        v_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
        # Surface Pressures
        p_surf = zeros(TF, nbodies, 2 * npanels - 1, naoa)
    end

    smooth_nodes = zeros(TF, nbodies, 2 * npanels - 1, naoa)

    ##### ----- Loop Through Bodies ----- #####
    for m in 1:nbodies

        # Panel Values
        # length
        panel_length = panels[m].panel_length
        # vector
        panel_vector = panels[m].panel_vector
        # edge locations
        panel_edges = panels[m].panel_edges
        # midpoints
        panel_midpoints = (panel_edges[:, 2, :] .+ panel_edges[:, 1, :]) ./ panel_length

        for a in 1:naoa

            ### --- Get Surface Distributions --- ###
            # - Get Surface Velocity - #
            #= NOTE:
                For xfoil-like method, the "solution" values ARE the velocity components on the surface.
            =#
            vti = [
                gamma0[i] * cosd(flow_angle[a]) + gamma90[i] * sind(flow_angle[a]) for
                i in nidx[m]
            ]

            # - Get Surface Pressure (Steady State) - #
            cpi = 1.0 .- vti .^ 2

            # Get Mean Pressure at PANEL MIDPOINTS
            cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0

            # organize outputs based on debug flag
            if debug
                v_surf[m][:, a] = vti
                p_surf[m][:, a] = cpi
                # Get Mean Pressure at PANEL MIDPOINTS
                cpibar = (cpi[1:(end - 1)] .+ cpi[2:end]) ./ 2.0
            else
                #smooth_distributions functions are found in utils.jl
                v_surf[m, :, a], smooth_nodes[m, :, a] = smooth_distributions(
                    Linear(), panel_edges, vti, npanels
                )
                p_surf[m, :, a], _ = smooth_distributions(
                    Linear(), panel_edges, cpi, npanels
                )
            end

            #---------------------------------#
            #      Calculate Coefficients     #
            #---------------------------------#

            #quarter chord location (moment reference location for inviscid case)
            x0 = chord / 4.0
            z0 = 0.0 #x0*sind(flow_angle[a]) #TODO should this be zero, or rotated with the airfoil?

            panelidx = mesh.mesh2panel

            ### --- Calculate Lift Coefficient --- ###
            cl[m, a] =
                sum([
                    cpibar[i] * (
                        -sind(flow_angle[a]) * panel_vector[i, 2] -
                        cosd(flow_angle[a]) * panel_vector[i, 1]
                    ) for i in panelidx[pidx[m]]
                ]) / chord

            ### --- Calculate Drag Coefficients --- ###
            #= NOTE:
                For the inviscid case, cdp is zero and cdi=cd for now.
            =#
            cdp[m, a] = 0.0
            #TODO: check of drag calculation needs to be negative or not
            cdi[m, a] =
                -sum([
                    cpibar[i] * (
                        cosd(flow_angle[a]) * panel_vector[i, 2] -
                        sind(flow_angle[a]) * panel_vector[i, 1]
                    ) for i in panelidx[pidx[m]]
                ]) / chord
            cd[m, a] = cdi[m, a]

            ### --- Calculate Moment Coefficient --- ###
            # initialize pieces of moment calculation
            cmmat = [2 1; 1 2] ./ 6.0

            # Moment arms
            dxddmi =
                panel_vector[panelidx[pidx[m]], 1] .*
                (panel_edges[panelidx[pidx[m]], 1, 1] .- x0) .+
                panel_vector[panelidx[pidx[m]], 2] .*
                (panel_edges[panelidx[pidx[m]], 1, 2] .- z0)

            dxddmip1 =
                panel_vector[panelidx[pidx[m]], 1] .*
                (panel_edges[panelidx[pidx[m]], 2, 1] .- x0) .+
                panel_vector[panelidx[pidx[m]], 2] .*
                (panel_edges[panelidx[pidx[m]], 2, 2] .- z0)

            # Moment Coefficient Calculation
            cm[m, a] =
                sum([
                    ([cpi[i] cpi[i + 1]] * cmmat * [
                        dxddmi[i]
                        dxddmip1[i]
                    ])[1] for i in panelidx[pidx[m]]
                ]) / chord^2
        end
    end

    # Again, returns depend on debug flag
    if debug
        return cl, cd, cdp, cdi, cm, v_surf, p_surf
    else
        return PlanarPolar(cl, cd, cdp, cdi, cm, v_surf, p_surf, smooth_nodes)
    end
end

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
