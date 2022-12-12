using FLOWFoil
using PyPlot

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

# c
c1 = a1 / k1
c2 = a2 / k2

# wedge angle (radians)
tau1 = 0.17279
tau2 = 0.17279

# get circle centers relative to their respective axes
center1 = [a1 * sin(beta1); a1 * cos(beta1) - c1]
center2 = [a2 * sin(beta2); a2 * cos(beta2) - c2]

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

figure(1)
clf()
plot(x1, y1)
plot([-1.5; 1.5], [0; 0], "k")
plot([0; 0], [-1.5; 1.5], "k")
plot([center1[1]], [center1[2]], "xC0"; markersize=7)

plot(x2, y2)
plot([1.9 - 0.75; 1.9 + 0.75], [-0.4; -0.4], "k")
plot([1.9; 1.9], [-0.4 - 0.75; -0.4 + 0.75], "k")
plot([center2[1] + 1.9], [center2[2] - 0.4], "xC1"; markersize=7)
axis("equal")

savefig("circles.pdf"; bbox_inches="tight")

################
### AIRFOILS ###
################

x1, z1 = FLOWFoil.karman_trefftz(
    center1, tau1 * 180.0 / pi; radius=k2, beta=beta1, normalize=false
)
x2, z2 = FLOWFoil.karman_trefftz(
    center2, tau2 * 180.0 / pi; radius=k2, beta=beta2, normalize=false
)
x2 .*= a2
z2 .*= a2
x2 .+= axis2[1]
z2 .+= axis2[2]

figure(2)
clf()
plot(x1, z1; label="circle 1")
plot(x2, z2; label="circle 2")
axis("equal")
