#=

It looks like there's a bug in the karman_trefftz parameterization that needs to be addressed.

Once fixed, this file will generate the geometry needed for the inviscid multi-airfoil validation for the mfoil implementation

=#

using FLOWFoil
using Plots

# axes to which circle geometry is relative
axis1 = [0.0; 0.0]
axis2 = [1.9; -0.4]

# circle radii
a1 = 1.096
a2 = 0.5

# angle between radius on x axis
beta1 = 0.001
beta2 = 0.08725

# k = a/c, a = radius
k1 = 1.096
k2 = 1.15

# wedge angle (radians)
tau1 = 0.17279
tau2 = 0.17279

###############
### CIRCLES ###
###############

N = 360
function getcircle(center, R, beta, N)
    t = range(-beta; stop=2 * pi - beta, length=N)
    x = R * cos.(t) .+ center[1]
    y = R * sin.(t) .+ center[2]
    return x, y
end

x1, y1 = getcircle(center1, a1, beta1, N)
x2, y2 = getcircle(center2, a2, beta2, N)
x2 .+= axis2[1]
y2 .+= axis2[2]

plot(x1, y1; aspectratio=1)
plot!([-1.5; 1.5], [0; 0]; color=:black)
plot!([0; 0], [-1.5; 1.5]; color=:black)
plot!([center1[1]], [center1[2]]; seriestype=:scatter, markercolor=1, markersize=7)

plot!(x2, y2; color=2)
plot!([1.9 - 0.75; 1.9 + 0.75], [-0.4; -0.4]; color=:black)
plot!([1.9; 1.9], [-0.4 - 0.75; -0.4 + 0.75]; color=:black)
plot!(
    [center2[1] + 1.9], [center2[2] - 0.4]; seriestype=:scatter, markercolor=2, markersize=7
)

################
### AIRFOILS ###
################

x1, z1 = FLOWFoil.AirfoilTools.karman_trefftz(beta1, k1, tau1; normalize=false)
x1 .*= a1
x2 .*= a1

x2, z2 = FLOWFoil.AirfoilTools.karman_trefftz(beta2, k2, tau2; normalize=false)
coords2 = [x2 z2]
AirfoilTools.rotate_coordinates!(coords2, 30.0; rotation_point=[0.0; 0.0])
coords2 .*= a2
coords2[:, 1] .+= axis2[1]
coords2[:, 2] .+= axis2[2]

plot(x1, z1; label="circle 1", aspectratio=1)
plot!(coords2[:, 1], coords2[:, 2]; label="circle 2")
