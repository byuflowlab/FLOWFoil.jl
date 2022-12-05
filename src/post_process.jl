#=
General Post Processing Types and Functions

Authors: Judd Mehr,
=#

abstract type Polar end

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
 - `lift::Float` : Lift Coefficient.
 - `drag::Float` : Total Drag Coefficient.
 - `pdrag::Float` : Pressure Drag Coefficient.
 - `idrag::Float` : Induced Drag Coefficient.
 - `moment::Float` : Moment Coefficient.
 - `surfacevelocity::Vector{Float}` : surface velocity distribution
 - `surfacepressure::Vector{Float}` : surface pressure distribution
"""
struct PlanarPolar{TF} <: Polar
    lift::Vector{TF}
    drag::Vector{TF}
    pdrag::Vector{TF}
    idrag::Vector{TF}
    moment::Vector{TF}
    surface_velocity::Array{TF,3}
    surface_pressure::Array{TF,3}
    xsmooth::Array{TF,3}
end

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#

"""

if debug = true, the outputs will not be in a struct format to allow for potentially non-equal lengths in the surface velocity and pressures
"""
function post_process(pp::PlanarProblem, problem, panels, mesh, solution; debug=false)
    return post_process(pp.singularity, problem, panels, mesh, solution; debug=debug)
end

"""
"""
function post_process(v::Vortex, problem, panels, mesh, solution; debug=false)
    return post_process_vortex(v.order, problem, panels, mesh, solution; debug=false)
end

"""
"""
function post_process_vortex(
    ::Linear, problem, panels, mesh, solution; npanels=80, debug=false
)

    ##### ----- Set Up ----- #####

    ### --- Rename for Convenience --- ###

    # number of airfoils
    nbodies = mesh.nbodies

    # number of angles of attack
    alpha = problem.angle_of_attack
    naoa = length(alpha)

    # node indices
    pidx = mesh.panel_indices
    nidx = mesh.node_indices

    # vortex strengths
    gamma0 = solution.x[:, 1]
    gamma90 = solution.x[:, 2]

    # chord length
    chord = mesh.chord

    # Panel Values
    # length
    panel_length = panels.panel_length
    # vector
    panel_vector = panels.panel_vector
    # edge locations
    panel_edges = panels.panel_edges
    # midpoints
    panel_midpoints = (panel_edges[:, 2, :] .+ panel_edges[:, 1, :]) ./ panel_length

    ### --- Initialize Outputs --- ###

    # output floating point type
    TF = typeof(chord)

    # - Coefficients - #

    # Lift coefficient
    cl = zeros(TF, nbodies)

    # Total drag coefficient
    cd = zeros(TF, nbodies)
    # Profile drag coefficient
    cdp = zeros(TF, nbodies)
    # Inviscid drag coefficient
    cdi = zeros(TF, nbodies)

    # Moment coefficient
    cm = zeros(TF, nbodies)

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
        for a in 1:naoa

            ### --- Get Surface Distributions --- ###
            # - Get Surface Velocity - #
            #= NOTE:
                For xfoil-like method, the "solution" values ARE the velocity components on the surface.
            =#
            vti = [
                gamma0[i] * cosd(alpha[a]) + gamma90[i] * sind(alpha[a]) for i in nidx[m]
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
            z0 = 0.0 #x0*sind(alpha[a]) #TODO should this be zero, or rotated with the airfoil?

            ### --- Calculate Lift Coefficient --- ###
            cl[m, a] =
                sum([
                    cpibar[i] * (
                        -sind(alpha[a]) * panel_vector[i, 2] -
                        cosd(alpha[a]) * panel_vector[i, 1]
                    ) for i in pidx[m]
                ]) / chord

            ### --- Calculate Drag Coefficients --- ###
            #= NOTE:
                For the inviscid case, cdp is zero and cdi=cd for now.
            =#
            cdp[m, a] = 0.0
            cdi[m, a] =
                sum([
                    cpibar[i] * (
                        cosd(alpha[a]) * panel_vector[i, 2] -
                        sind(alpha[a]) * panel_vector[i, 1]
                    ) for i in pidx[m]
                ]) / chord
            cd[m, a] = cdi[m, a]

            ### --- Calculate Moment Coefficient --- ###
            # initialize pieces of moment calculation
            cmmat = [2 1; 1 2] ./ 6.0

            # Moment arms
            dxddmi =
                panel_vector[pidx[m], 1] .* (panel_edges[pidx[m], 1, 1] .- x0) .+
                panel_vector[pidx[m], 2] .* (panel_edges[pidx[m], 1, 2] .- z0)

            dxddmip1 =
                panel_vector[pidx[m], 1] .* (panel_edges[pidx[m], 2, 1] .- x0) .+
                panel_vector[pidx[m], 2] .* (panel_edges[pidx[m], 2, 2] .- z0)

            # Moment Coefficient Calculation
            cm[m, a] =
                sum([
                    ([cpi[i] cpi[i + 1]] * cmmat * [dxddmi[i]; dxddmip1[i]])[1] for
                    i in pidx[m]
                ]) / chord^2
        end
    end

    # Again, returns depend on debug flag
    if debug
        return sum(cl), sum(cd), sum(cdp), sum(cdi), sum(cm), v_surf, p_surf
    else
        return PlanarPolar(cl, cd, cdp, cdi, cm, v_surf, p_surf, smooth_nodes)
    end
end

"""
"""
function smooth_distributions(::Linear, panel_edges, distribution, npanels)

    #= NOTE:
        Akima splines in FLOWMath require the 'x' values to be monotonically ascending.
        Therefore, we need to get all the panel edge points and then divide them into top and bottom in order to create our splines.
    =#

    # - Get 'x' values from panel edges - #
    x = [panel_edges[:, 1, 1]; panel_edges[end, 2, 1]]

    # - Split the 'x' values - #

    # find the minimum and index
    minx, minidx = findmin(x)

    # the bottom needs to be flipped to ascend monotonically
    xbot = x[minidx:-1:1]
    # the top is already in the right direction
    xtop = x[minidx:end]

    # - Get smooth 'x' values from cosine spacing - #
    # Get cosine spaced values from zero to one.
    xcosine = cosine_spacing(npanels)

    # - Transform the cosine spaced values to the minimum and maximum points - #
    # Get the maximum x value
    maxx = maximum(x)

    xsmooth = linear_transform([0.0; 1.0], [minx; maxx], xcosine)

    # - Generate smooth distribution - #
    distbot = FLOWMath.akima(xbot, distribution[minidx:-1:1], xsmooth)
    disttop = FLOWMath.akima(xtop, distribution[minidx:end], xsmooth)

    # - Combine distribution and x values - #
    xs = [reverse(xsmooth); xsmooth[2:end]]
    dist = [reverse(distbot); disttop[2:end]]

    # - Return - #
    return dist, xs
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

#---------------------------------#
#            FUNCTIONS            #
#---------------------------------#
