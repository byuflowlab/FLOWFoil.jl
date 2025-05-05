using Plots, Revise, FLOWFoil, .AirfoilTools, FLOWMath

function close_te(x, y; n_smoothing=5)
    x_closed = copy(x)
    y_closed = copy(y)

    smooth_factor = FLOWMath.linear(
        x[(end - n_smoothing):end], range(1, 0, n_smoothing + 1), x[(end - n_smoothing):end]
    )

    for (i, sf) in enumerate(smooth_factor)
        avg_te_x = (x[i] + x[end - i + 1]) / 2.0
        avg_te_y = (y[i] + y[end - i + 1]) / 2.0
        x_closed[i] += abs(x_closed[i] - avg_te_x) * sf
        x_closed[end - i + 1] -= abs(x_closed[end - i + 1] - avg_te_x) * sf
        y_closed[i] += abs(y_closed[i] - avg_te_y) * sf
        y_closed[end - i + 1] -= abs(y_closed[end - i + 1] - avg_te_y) * sf
    end

    return x_closed, y_closed
end

#x,y = naca65_scaled(1.2)
#x, y = close_te(x, y)
naca4_parameters = NACA4(2.0, 5.0, 12.0, false)
naca4_parameters2 = NACA4(0.0, 0.0, 16.0, false)
naca4_parameters3 = NACA4(1.0, 4.0, 12.0, false)
num_points = 161
x,y = naca4(naca4_parameters, N = num_points)
x2, y2 = naca4(naca4_parameters2, N = num_points)
x3,y3 = naca4(naca4_parameters3, N = num_points)

#linear transform for body 1 ### only for Lewis case
#y = y .+ 1.0

#linear transform for body 2
theta = 20
x2 = x2.*cosd(theta) .+ y2.*sind(theta)
y2 = -x2.*sind(theta) .+ y2.*cosd(theta)
x2 = x2.*0.5
y2 = y2.*0.5
x2 = x2 .+ 1.0
y2 = y2 .- 0.05
#y2 = y2 .+ 1.0 ##only use this line for lewis

#linear transform for body 3
theta2 = 20
x3 = x3.*cosd(theta2) .+ y3.*sind(theta2)
y3 = -x3.*sind(theta2) .+ y3.*cosd(theta2)
x3 = x3.*0.25
y3 = y3.*0.25
x3 = x3 .+ 1.5
y3 = y3 .- 0.25

#compile vectors into tuple
a = [x y]
b = [x2 y2]
c = [x3 y3]
coordinates = (a, b)


flow_angles = [-1.0, 0.0, 1.0]
flow_angles_lewis = [0.0]
method1 = Mfoil(false)
method2 = FLOWFoil.Martensen(false, 0.001, 0.0, 10.0, false)
method3 = FLOWFoil.Lewis([false])
method4 = FLOWFoil.HessSmith(1.0)
sol = FLOWFoil.analyze(x,y, flow_angles, method = method4)
println(sol.cl)


##### Plots ###########################

#plot multi body case
#=
plot(x,y, aspectratio = 1, label = "Body 1", grid = false)
plot!(x2, y2, label = "Body 2")
plot!(x3, y3, label = "Body 3")
=#


##Martensen Grid Convergence Study - multibody
#=
cl_body_1 = [0.5078, 0.3171, 0.2936, 0.2686, 0.2615]
cl_body_2 = [-0.0424,-0.0264, -0.02396, -00.0213, -0.02057]
num_points_vector = [161, 362, 500, 1000, 1500]
plot(num_points_vector, cl_body_1, grid = false, label = "Body 1", xlabel = "Number of Points", ylabel = "Cl")
plot!(num_points_vector, cl_body_2, label = "Body 2")
=#

#MFoil Grid Convergence Study - multibody
#=
cl_body_1 = [-2.0704, -2.1634, -10.4414, -4.1817, -1.8214]
cl_body_2 = [3.8753, 0.9787, 116.454, -1.1943, -1.5972]
num_points_vector = [161, 362, 500, 1000, 1500]
plot(num_points_vector, cl_body_1, grid = false, label = "Body 1", xlabel = "Number of Points", ylabel = "Cl", ylims=(-5.0,5.0))
plot!(num_points_vector, cl_body_2, label = "Body 2")
=#
###################################################3
