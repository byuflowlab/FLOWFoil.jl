function generate_system_matrices(method::Martensen, panels::TP, mesh) where {TP<:Panel}
    return generate_system_matrices(method, [panels], mesh)
end

function generate_system_matrices(method::Martensen, panels, mesh)
    # Get coeffiecient matrix (A, left hand side)
    A = assemble_coupling_matrix(panels, mesh, method.pitch, method.stagger)
    #assemble_periodic_influence_matrix(method.singularity, panels, mesh) #not sure how this function works so I commented it out

    # Get boundary conditions (b, right hand side)
    b = assemble_periodic_boundary_conditions(method.boundary, panels, mesh)

    return InviscidSystem(A, b, mesh.panel_indices)
end

#---------------------------------#
#       COEFFICIENT MATRIX        #
#---------------------------------#

#this function essentially takes in the panel geometry and outputs the coupling matrix
#Important: panels need to include cascade pitch 
function assemble_coupling_matrix(panels::TF, mesh, pitch, stagger)
    m = panels.npanels
    coup = Array{TF, 2}(undef, m, m) .= 0.0 #coup defined using number of panels from panel_geometry.jl
    
    #compute self-inducing coupling coefficients
    coup[1,1] = -0.5 - (panels.panel_angle[2] - panels.panel_angle[m] - 2*pi) / (8*pi)
    coup[m,m] = -0.5 - (panels.panel_angle[1] - panels.panel_angle[m - 1] - 2*pi) / (8*pi)
    for i = 2:m - 1
        coup[i,i] = -0.5 - (panels.panel_angle[i + 1] - panels.panel_angle[i - 1]) / (8*pi)
    end

    #add stagger to the cosine and sine vectors to follow Lewis
    for i = 1:m
        #here I re-evaluate the sine and cosine vectors
        panels.sine_vector[i] = sin(stagger + panels.panel_angle[i])
        panels.cosine_vector[i] = cos(stagger + panels.panel_angle[i])
    end
    #define convenience variables
    xm = 0.0
    ym = 0.0
    xn = 0.0
    yn = 0.0
    chord = 1.0

    for i = 1:m
        for j = i:m
            if j != i
                #if pitch / chord is greater than 30, then the program will perform a single airfoil analysis
                if pitch / chord > 30.0
                    xm = panels.panel_center[i]
                    xn = panels.panel_center[j]
                    ym = panels.panel_center[i]
                    yn = panels.panel_center[j]
                    r = (xn - xm)^2 + (yn -m)^2
                    u = (yn - ym) / (r*2*pi)
                    v = -(xn - xm) / (r*2*pi)
                    coup[j,i] = (u*panels.cosine_vector[j] + v*panels.sine_vector[j])*panels.panel_length[i]
                    coup[i,j] = -(u*cosine_vector[i] + v*sine_vector[i])*panels.panel_length[j]
                else
                    xm = panels.panel_center[i]
                    xn = panels.panel_center[j]
                    ym = panels.panel_center[i]
                    yn = panels.panel_center[j]
                    a = ((xm - xn)*cos(stagger) - (ym - yn)*sin(stagger))*2*pi / pitch
                    b = ((xm - xn)*sin(stagger) + (ym - yn)*cos(stagger))*2*pi / pitch
                    k = 0.5 / pitch / (cosh(a) - cos(b))
                    coup[j,i] = (sinh(a)*panels.sine_vector[j] - sin(b)*panels.cosine_vector[j])*k*panels.panel_length[i]
                    coup[i,j] = (-sinh(a)*sine_vector[i] + sin(b)*cosine_vector[i])*k*panels.panel_length[j]
                end
            end
        end
    end

    #apply the back-diagonal correction
    sum = 0.0
    for i = 1:m
        sum = 0.0
        for j = 1:m
            if j != m + 1 - i
                sum = sum + coup[j, i]*panels.panel_length[j]
            end
            coup[m + 1 - i, i] = -sum / panels.panel_length[m + 1 - i]
        end
    end
end

"""
    assemble_periodic_influence_matrix(v::Vortex, body_of_revolution, panels, mesh)

Assembles the "A" matrix (left hand side coefficient matrix).

# Arguments:
- `s::Singularity` : The singularity type used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

# Returns:
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_periodic_influence_matrix(::Singularity, panels, mesh) end

function assemble_periodic_influence_matrix(v::Vortex, panels, mesh)
    return assemble_periodic_vortex_matrix(v.order, panels, mesh)
end

"""
    assemble_periodic_vortex_matrix(::Order,  panels, mesh)

Assembles the coefficient matrix for a given order of singularity.

# Arguments:
- `o::Order` : The order of singularity used.
- `:Vector{Bool}` : flags whether bodies are bodies of revolution or not.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

# Returns:
- A::Matrix{Float}` : The influence coefficient matrix for the linear system
"""
function assemble_periodic_vortex_matrix(::Order, panels, mesh) end

function assemble_periodic_vortex_matrix(::Constant, panels, mesh)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies

    # initialize coefficient matrix
    TF = eltype(mesh.x)
    amat = zeros(TF, (N + nbodies, N + nbodies))

    # Loop through system

    ### --- Loop through bodies --- ###
    for m in 1:nbodies
        for n in 1:nbodies
            ### --- Loop through panels --- ###
            for i in idx[m]
                for j in idx[n]

                    ### --- Calculate influence coefficient --- ###
                    s = j >= idx[m][end] - i ? -1.0 : 1.0
                    amat[i, j] =
                        s * calculate_periodic_vortex_influence(
                            Constant(), panels[m], panels[n], mesh, i, j
                        )
                end
            end

            if m == n

                ### --- Apply Back Substitution --- ###
                for i in idx[n]
                    sum = 0.0
                    jidx = idx[n][end] + 1 - i
                    for j in idx[m]
                        if j != jidx
                            sum += amat[j, i] * panels[n].panel_length[j]
                        end
                    end
                    dmagj = panels[n].panel_length[jidx]
                    amat[jidx, i] = -sum / dmagj
                end

                ### --- Apply Kutta Condition --- ###
                # put in the kutta condition for each airfoil (end rows of the system matrix)
                amat[N + m, idx[n][1]] = 1.0
                amat[N + m, idx[n][end]] = 1.0

                #put unit bound vortex value in each row
                amat[idx[m], idx[n][end] + n] .= 1.0
            end

            # # NOTE: this doesn't seem to change anything...
            # # Bound Vortex Correction
            # for i in 1:N
            #     for j in 1:M
            #         amat[i, j] += meshj.panels[j].length
            #     end
            # end

        end
    end

    return amat
end

#---------------------------------#
#    BOUNDARY CONDITION MATRIX    #
#---------------------------------#

"""
    assemble_periodic_boundary_conditions(::BoundaryCondition,  panels, mesh)

Assemble boundary condition vector.

# Arguments:
- `bc::BoundaryCondition` : The type of boundary condition to be used.
- `panels::Vector{Panel}` : Vector of panel objects (one for each body in the system)
- `mesh::Mesh` : The mesh object containing relative geometry for the influence coefficient calculations.

# Returns
 - `b::Matrix{Float}` : Boundary condition matrix
"""
function assemble_periodic_boundary_conditions(::BoundaryCondition, panels, mesh) end

function assemble_periodic_boundary_conditions(::Neumann, panels, mesh)

    ### --- SETUP --- ###

    # - Rename for Convenience - #
    idx = mesh.panel_indices
    N = idx[end][end]
    nbodies = mesh.nbodies

    # initialize boundary condition array
    TF = eltype(mesh.x)
    bc = zeros(TF, N + nbodies, 3)

    ### --- Loop through bodies --- ###
    for m in 1:nbodies

        # generate portion of boundary condition array associated with mth body
        bc[idx[m], 1] = [-cos(panels[m].panel_angle[i]) for i in idx[m]]
        bc[idx[m], 2] = [-sin(panels[m].panel_angle[i]) for i in idx[m]]
        bc[idx[m], 3] = [1.0 for i in idx[m]]
    end

    return bc
end

######################################################################
#                                                                    #
#                        FROM CHATGPT TRANSLATION OF LEWIS                         #
#                                                                    #
######################################################################

function coupling_coefficients()
    # Variable declarations
    r = 0.0
    u = 0.0
    v = 0.0
    a = 0.0
    b = 0.0
    k = 0.0
    sinh = 0.0
    cosh = 0.0
    e = 0.0
    pitch = 0.0
    chord = 0.0
    m = 0
    stagger = 0.0
    slope = zeros(Float64, m)
    xdata = zeros(Float64, m)
    ydata = zeros(Float64, m)
    ds = zeros(Float64, m)
    coup = zeros(Float64, m, m)

    # Arrays to store sine and cosine values
    sine = zeros(Float64, m)
    cosine = zeros(Float64, m)

    # Loop to initialize sine and cosine arrays
    for i in 1:m
        sine[i] = sin(stagger + slope[i])
        cosine[i] = cos(stagger + slope[i])
    end

    # Self inducing coupling coefficients
    coup[1, 1] = -0.5 * (slope[2] - slope[m] - 2.0 * pi) / (8.0 * pi)
    coup[m, m] = -0.5 * (slope[1] - slope[m - 1] - 2.0 * pi) / (8.0 * pi)

    for i in 2:(m - 1)
        coup[i, i] = -0.5 * (slope[i + 1] - slope[i - 1]) / (8.0 * pi)
    end

    # Nested loops for coupling coefficients
    for i in 1:m
        for j in 1:m
            if j < i
                if pitch / chord > 30.0
                    # Revert to single aerofoil for very wide blade spacing
                    r = sqrt((xdata[i] - xdata[j])^2 + (ydata[i] - ydata[j])^2)
                    u = (ydata[j] - ydata[i]) / (r * 2 * pi)
                    v = -(xdata[j] - xdata[i]) / (r * 2 * pi)
                    coup[i, j] = (u * cosine[j] + v * sine[j]) * ds[i]
                    coup[j, i] = -(u * cosine[i] + v * sine[i]) * ds[j]
                else
                    # Cascade coupling coefficients
                    a =
                        (
                            (xdata[i] - xdata[j]) * cos(stagger) -
                            (ydata[i] - ydata[j]) * sin(stagger)
                        ) *
                        2 *
                        pi / pitch
                    b =
                        (
                            (xdata[i] - xdata[j]) * sin(stagger) +
                            (ydata[i] - ydata[j]) * cos(stagger)
                        ) *
                        2 *
                        pi / pitch
                    e = exp(a)
                    sinh = 0.5 * (e - 1.0 / e)
                    cosh = 0.5 * (e + 1.0 / e)
                    k = 0.5 * pitch / (cosh - cos(b))
                    coup[i, j] = (sinh * sine[i] - sin(b) * cosine[i]) * k * ds[i]
                    coup[j, i] = -(sinh * sine[j] + sin(b) * cosine[j]) * k * ds[j]
                end
            end
        end
    end
end

