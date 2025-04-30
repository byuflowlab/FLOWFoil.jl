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

#=
x,y = naca65_scaled(1.2)
x, y = close_te(x, y)
=#
naca4_parameters = NACA4(2.0, 4.0, 12.0, false)
x,y = naca4(naca4_parameters)
flow_angles = [-1.0, 0.0, 1.0]
reynolds = [1e6]
machs = [0.0]
method1 = Mfoil()
sol = FLOWFoil.analyze(x,y, flow_angles, reynolds, machs, method = method1)