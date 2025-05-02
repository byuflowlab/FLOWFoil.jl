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
x,y = naca4(naca4_parameters)
x2, y2 = naca4(naca4_parameters2)
theta = 30
x2 = x2.*cosd(theta) .+ y2.*sind(theta)
y2 = -x2.*sind(theta) .+ y2.*cosd(theta)
x2 = x2.*0.5
y2 = y2.*0.5
x2 = x2 .+ 1.0
y2 = y2 .- 0.05
a = [x y]
b = [x2 y2]
coordinates = (a, b)
flow_angles = [-1.0, 0.0, 1.0]
method1 = Mfoil(false)
method2 = FLOWFoil.Martensen(true, 1.0, 0.0, 10.0, false)
method3 = FLOWFoil.Lewis([false])
sol = FLOWFoil.analyze(coordinates, flow_angles)
println(sol.cl[:,2])
plot(x,y, aspectratio = 1)
plot!(x2,y2)
