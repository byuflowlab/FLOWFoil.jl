#TODO: update these functions and docstrings

"""
    rawkarmanTrefftzFoil(center::Array{Float64}=[-0.1 0.1], wedge::Float64=0.0;N::Int=360)

Calculate geometry of Karman-Trefftz airfoil based on the center point of a cylinder and a wedge angle.
"""
function rawkarmantrefftz(center=[-0.1 0.1], wedge=0.0; N=360)
    r = sqrt((1 - center[1])^2 + center[2]^2) #radius
    beta = asin(center[2] / r) #solve for beta
    #     #Stagnation Points
    #     A = center[1]-r*cos(beta)
    #     B = center[1]+r*cos(beta)
    wedge *= pi / 180 #convert to radians
    lambda = 2 - wedge / pi

    #initialize
    theta = range(0.0; stop=2 * pi, length=N) #some say this needs to go from (-beta, 2*pi-beta)
    w = zeros(N) * im

    #solve for geometry
    zeta(theta) = 1 .+ r * (exp(im * theta) - exp(-im * beta))
    inp(theta) = (1 .+ 1 / zeta(theta))^(lambda)
    inm(theta) = (1 .- 1 / zeta(theta))^(lambda)
    for i in 1:N
        w[i] = lambda * (inp(theta[i]) + inm(theta[i])) / (inp(theta[i]) - inm(theta[i])) #karman-trefftz transformation
    end

    #separate x and z coordinates
    x = real(w)
    z = imag(w)

    return x, z, theta
end

"""
    karmantrefftz(center::Array{Float64}=[-0.1 0.1], wedge::Float64=0.0;N::Int=360)

Calculate geometry of Karman-Trefftz airfoil based on the center point of a cylinder and a wedge angle.
"""
function karmantrefftz(center=[-0.1 0.1], wedge=0.0; N=360)
    x, z = rawkarmantrefftz(center, wedge; N=N) #get coordinates
    minx = minimum(x) #find minimum x (probably around -2)
    maxx = maximum(x) #find maximum x (probably around +2)
    chord = maxx - minx #find chord
    x .+= abs(minx) #translate airfoil so that leading edge is at zero
    x ./= chord #normalize x coordinates by chord length
    z ./= chord #normalize y coordinates by chord length
    #cleanup
    oneloc = findall(x -> x == 1.0, x) #find first 1.0 in x array
    if oneloc[1] > 1
        shiftval = length(x) - oneloc[1] + 1
        x = circshift(x, shiftval)
        z = circshift(z, shiftval)
        theta = circshift(theta, shiftval)
    end

    #Split
    LEidx = argmin(x)
    xu = x[LEidx:-1:1]
    zu = z[LEidx:-1:1]
    thetau = theta[LEidx:-1:1]
    xl = x[(LEidx + 1):1:end]
    zl = z[(LEidx + 1):1:end]
    thetal = theta[(LEidx + 1):1:end]

    return xu, zu, thetau, xl, zl, thetal, chord
end

"""
    rawjoukowsky(center::Array{Float64}=[-0.1 0.1] ; N::Int=360)

Calculate geometry of Joukowsky Airfoil based on a cylinder's center point.
"""
function rawjoukowsky(center::Array{Float64}=[-0.1 0.1]; N::Int=360)
    r = sqrt((1 - center[1])^2 + center[2]^2) #radius #rename for convenience
    beta = asin(center[2] / r) #solve for beta
    #     #Stagnation Points
    #     A = center[1]-r*cos(beta)
    #     B = center[1]+r*cos(beta)
    #     println("Stagnations: \n", A,"\n",B)
    #initialize
    theta = range(0.0; stop=2 * pi, length=N) #some say this needs to go from (-beta, 2*pi-beta)
    w = zeros(N) * im

    #Solve for geometry
    zeta(theta) = 1 .+ r * (exp(im * theta) - exp(-im * beta)) #WHY??? Where in the world does this come from?
    for i in 1:360
        w[i] = zeta(theta[i]) + 1 / zeta(theta[i]) #Joukowsky transformation w = z +1/z
    end

    #separate x and z coordinates
    x = real(w)
    z = imag(w)

    return x, z, theta, r, beta
end

function rawjoukowsky2(center, R; N=360)
    beta = asin(center[2] / R) #solve for beta

    theta = reverse(range(-beta; stop=2 * pi - beta, length=N))

    a = center[1] + R * cos(beta)

    xi =
        (1.0 .+ R / a * (exp.(im * theta) .- exp(-im * beta))) .+
        1.0 ./ (1.0 .+ R / a * (exp.(im * theta) .- exp(-im * beta)))
    xi *= a

    x = real(xi)
    y = imag(xi)

    return x, y, theta, beta, a
end

function joukowskysurface(center, R, alpha, U; N=360)

    x, y, theta, beta, a = rawjoukowsky2(center, R; N=N)

    z = a * (1.0 .+ R / a * (exp.(im * theta) .- exp(-im * beta)))

    vmag =
        2.0 * U * (sin.(theta) .- alpha .+ sin(alpha) .+ beta) .* abs.(z) ./
        (abs.((z .- a^2) ./ z))

    cp = 1.0 .- vmag .^ 2.0 / U^2.0

    return x, y, vmag, cp
end

"""
    joukowsky(center::Array{Float64}=[-0.1 0.1] ; N::Int=360)

Calculate the unit Joukowsky Airfoil based on a cylinder's center point.
"""
function joukowsky(center::Array{Float64}=[-0.1 0.1], ; N::Int=360)
    x, z = rawjoukowsky(center; N=N) #get coordinates
    minx = minimum(x) #find minimum x (probably around -2)
    maxx = maximum(x) #find maximum x (probably around +2)
    chord = maxx - minx #find chord length
    x .+= abs(minx) #shift airfoil so that leading edge is at zero
    x ./= chord #normalize x coordinates by chord length
    z ./= chord #normalize y coordinates by chord length
    #cleanup
    oneloc = findall(x -> x == 1.0, x) #find first 1.0 in x array
    if oneloc[1] > 1
        shiftval = length(x) - oneloc[1] + 1
        x = circshift(x, shiftval)
        z = circshift(z, shiftval)
    end

    #Split
    LEidx = argmin(x)
    xu = x[LEidx:-1:1]
    zu = z[LEidx:-1:1]
    xl = x[(LEidx + 1):1:end]
    zl = z[(LEidx + 1):1:end]

    return xu, zu, xl, zl
end
