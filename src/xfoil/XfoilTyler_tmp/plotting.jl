clr1 = (0.227,0.365,0.490)   # blue
clr2 = (0.478,0.557,0.329)   # green
clr3 = (0.549,0.376,0.247)   # brown
clr4 = (0.490,0.345,0.427)   # purple
clr5 = (0,0,0.)              # black
clr6 = (167, 9, 9) ./ 255    # red
clr7 = (212, 113, 47) ./ 255 # orange

function plotting_functions_sweep(alphas, panels_airfoil, inviscidoutputs, 
                            viscousoutputs, 
                            which_plots;
                            data_file=nothing)

    if :geometry in which_plots
        xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[end].x2)
        ys = vcat([p.y1 for p in panels_airfoil], panels_airfoil[end].y2)

        figure()
        ax = PyPlot.axes()
        plot(xs, ys, color=clr1, linewidth=2)
        axis("equal")
    end

    if :cl in which_plots
        if isnothing(viscousoutputs)
            cl = inviscidoutputs.cl
        else
            cl = viscousoutputs.cl
        end

        figure()
        plot(alphas*180/pi, cl, color=clr1, linewidth=2)
        xlabel("alpha (deg)")
        ylabel("cl")
    end

    if :cd in which_plots
        cd = viscousoutputs.cd

        figure()
        plot(alphas*180/pi, cd, color=clr1, linewidth=2)
        xlabel("alpha (deg)")
        ylabel("cd")
    end

    if :drag_polar in which_plots
        cl = viscousoutputs.cl
        cd = viscousoutputs.cd

        figure()
        plot(cd, cl, color=clr1, linewidth=2)
        xlabel("cd")
        ylabel("cl")
    end
end

function plotting_functions(panels_airfoil, inviscidoutputs, 
                            viscousoutputs, 
                            which_plots; 
                            data_file=nothing, verbose=true, 
                            alphas=nothing, Res=nothing, 
                            Machs=nothing, rhos=nothing,
                            Vinfs=nothing)

    close("all")

    if !isnothing(viscousoutputs)
        (; s, su, sl, sw, idx_u, idx_l, idx_w, stu, stl, s_stag, Vwake,
            deltas_u, thetas_u, state3s_u, 
            deltas_l, thetas_l, state3s_l, 
            deltas_w, thetas_w, state3s_w,
            idx_transition_l, idx_transition_u, cl, cd, 
            transition_upper, transition_lower,
            panels_wake, panel_TE) = viscousoutputs
    else
        (; cl, cd) = inviscidoutputs
    end

    if data_file != nothing
        include(joinpath("data/mfoil", data_file))
    
        println("")
        println("cl: \t\t$cl")
        if !isnothing(data_file)
            println("cl_mfoil: \t$cl_mfoil")
            try
                println("Relative error: $(Int(round(1000000*abs(cl_mfoil - cl)/cl_mfoil))/10000)%")
            catch
                println("Relative error: Too big!")
            end
        end
        println("")
        println("cd: \t\t$cd")
        if !isnothing(data_file)
            println("cd_mfoil: \t$cd_mfoil")
            try
                println("Relative error: $(Int(round(1000000*abs(cd_mfoil - cd)/cd_mfoil))/10000)%")
            catch
                println("Relative error: Too big!")
            end
        end
        println("")
        println("transition on upper surface: \t$transition_upper")
        println("transition on lower surface: \t$transition_lower")
        println("(measured from stagnation point)")
        println("")
    end

    # println("Transition point: $(viscousparameters.xt[1])")

    if :geometry in which_plots
        xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[end].x2)
        ys = vcat([p.y1 for p in panels_airfoil], panels_airfoil[end].y2)

        figure()
        ax = PyPlot.axes()
        if isnothing(viscousoutputs)
            plot(xs, ys, color=clr1, linewidth=2)
            # scatter(xs, ys, color=clr1, label="panel endpoints")

        else
            plot(xs[idx_u[1]:idx_u[end]], ys[idx_u[1]:idx_u[end]], color=clr1, linewidth=2)
            # scatter(xs[idx_u[1]:idx_u[end]], ys[idx_u[1]:idx_u[end]], color=clr1)
            plot(xs[idx_l[end]:idx_l[1]], ys[idx_l[end]:idx_l[1]], color=clr6, linewidth=2)
            # scatter(xs[idx_l[end]:idx_l[1]], ys[idx_l[end]:idx_l[1]], color=clr6)

            # stagnation
            idx_stag = FLOWMath.linear(s, 1:length(s), s_stag)
            idx_m1 = Int(floor(idx_stag))
            idx_p1 = Int(ceil(idx_stag))
            xs_stag = FLOWMath.linear([idx_m1, idx_p1], [xs[idx_m1], xs[idx_p1]], idx_stag)
            ys_stag = FLOWMath.linear([idx_m1, idx_p1], [ys[idx_m1], ys[idx_p1]], idx_stag)

            scatter([xs_stag, xs_stag], [ys_stag, ys_stag], color=clr4, marker="x", s=120, label="stagnation")

            #transition
            if stu <= su[end]
                idx_transition = FLOWMath.linear(su, idx_u[1]:idx_u[end], stu)
                idx_m1 = Int(floor(idx_transition))
                idx_p1 = Int(ceil(idx_transition))
                xs_transition = FLOWMath.linear([idx_m1, idx_p1], [xs[idx_m1], xs[idx_p1]], idx_transition)
                ys_transition = FLOWMath.linear([idx_m1, idx_p1], [ys[idx_m1], ys[idx_p1]], idx_transition)

                scatter([xs_transition, xs_transition], [ys_transition, ys_transition], color=clr2, marker="x", s=100, label="transition")
            else
                # plot at end
                scatter([xs[end], xs[end]], [ys[end], ys[end]], color=clr2, marker="x", s=100, label="transition")
            end

            if stl <= sl[end]
                idx_transition = FLOWMath.linear(sl, idx_l[1]:-1:idx_l[end], stl)
                idx_m1 = Int(floor(idx_transition))
                idx_p1 = Int(ceil(idx_transition))

                xs_transition = FLOWMath.linear([idx_m1, idx_p1], [xs[idx_m1], xs[idx_p1]], idx_transition)
                ys_transition = FLOWMath.linear([idx_m1, idx_p1], [ys[idx_m1], ys[idx_p1]], idx_transition)

                scatter([xs_transition, xs_transition], [ys_transition, ys_transition], color=clr2, marker="x", s=100)
            else
                # plot at end
                scatter([xs[1], xs[1]], [ys[1], ys[1]], color=clr2, marker="x", s=100)
            end

            if !isnothing(panels_wake)
                xs_wake = vcat([p.x1 for p in panels_wake], panels_wake[end].x2)
                ys_wake = vcat([p.y1 for p in panels_wake], panels_wake[end].y2)

                plot(xs_wake, ys_wake, color=clr5)
                # scatter(xs_wake, ys_wake, color=clr5)
            end
            legend()
        end

        if !isnothing(panel_TE)
            plot([panel_TE.x1, panel_TE.x2], [panel_TE.y1, panel_TE.y2], color=clr3, linewidth=2)
            # scatter([panel_TE.x1, panel_TE.x2], [panel_TE.y1, panel_TE.y2], color=clr3)
        end
        xlabel(L"x/c")
        ylabel(L"y/c")
        axis("equal")
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)
    end

    if :cp in which_plots
        xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[1].x1)
        # xs = panels.xbar 
        cps = isnothing(viscousoutputs) ? inviscidoutputs.cps : viscousoutputs.cps

        figure()
        ax = PyPlot.axes()
        if isnothing(viscousoutputs)
            plot(xs, -cps, color=clr1, linewidth=2, label="model")
        else
            plot(xs[idx_u[1]:idx_u[end]], -cps[idx_u[1]:idx_u[end]], color=clr1, linewidth=2, label="model")
            plot(xs[idx_l[1]:-1:idx_l[end]], -cps[idx_l[1]:-1:idx_l[end]], color=clr6, linewidth=2)
            
            if !isnothing(panels_wake)
                plot(xs_wake, -cps[idx_w[1]:idx_w[end]], color=clr5, linewidth=2)
            end
        end

        ylabel(L"-c_p")
        xlabel(L"x/c")
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        if data_file != nothing
            if isnothing(viscousoutputs)
                plot(x_mfoil, -cps_mfoil, "--", color=clr6, label="data")
            else
                plot(x_mfoil[idx_u_mfoil], -cps_mfoil[idx_u_mfoil], "--", color=clr1, label="data")
                plot(x_mfoil[idx_l_mfoil], -cps_mfoil[idx_l_mfoil], "--", color=clr6)
                plot(x_mfoil[idx_w_mfoil], -cps_mfoil[idx_w_mfoil], "--", color=clr5)
            end
            legend()
        end
    end

    if :cp_paper in which_plots
        xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[1].x1)
        # xs = panels.xbar 
        cps = isnothing(viscousoutputs) ? inviscidoutputs.cps : viscousoutputs.cps

        figure()
        ax = PyPlot.axes()

        if data_file != nothing
            plot(x_mfoil[idx_u_mfoil], -cps_mfoil[idx_u_mfoil], color=clr1, label="mfoil", linewidth=2)
            plot(x_mfoil[idx_l_mfoil], -cps_mfoil[idx_l_mfoil], color=clr1, linewidth=2)
            if !isnothing(viscousoutputs)
                plot(x_mfoil[idx_w_mfoil], -cps_mfoil[idx_w_mfoil], linewidth=2, color=clr1)
            end
            legend()
        end

        if isnothing(viscousoutputs)
            plot(xs, -cps, "--", color=clr6, linewidth=2, label="FLOWFoil")
        else
            plot(xs[idx_u[1]:idx_u[end]], -cps[idx_u[1]:idx_u[end]], color=clr6, "--", linewidth=2, label="FLOWFoil")
            plot(xs[idx_l[1]:-1:idx_l[end]], -cps[idx_l[1]:-1:idx_l[end]], color=clr6, "--", linewidth=2)
            
            if !isnothing(panels_wake)
                plot(xs_wake, -cps[idx_w[1]:idx_w[end]], color=clr6, "--", linewidth=2)
            end
        end

        ylabel(L"-c_p")
        xlabel(L"x/c")
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)
        legend()
    end

    if :cp_comparison in which_plots
        xs = vcat([p.x1 for p in panels_airfoil], panels_airfoil[1].x1)
        # xs = panels.xbar 
        cps_inviscid = inviscidoutputs.cps
        cps_viscous = viscousoutputs.cps

        figure()
        ax = PyPlot.axes()
        plot(xs[idx_u[1]:idx_u[end]], -cps_viscous[idx_u[1]:idx_u[end]], color=clr1, linewidth=2, label="viscous")
        plot(xs[idx_l[1]:-1:idx_l[end]], -cps_viscous[idx_l[1]:-1:idx_l[end]], color=clr6, linewidth=2)
        plot(xs[idx_u[1]:idx_u[end]], -cps_inviscid[idx_u[1]:idx_u[end]], "--", color=clr1, linewidth=2, label="inviscid")
        plot(xs[idx_l[1]:-1:idx_l[end]], -cps_inviscid[idx_l[1]:-1:idx_l[end]], "--", color=clr6, linewidth=2)
            
        if !isnothing(panels_wake)
            plot(xs_wake, -cps[idx_w[1]:idx_w[end]], color=clr5, linewidth=2)
        end

        ylabel(L"-c_p")
        xlabel(L"x/c")
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)
        legend()
    end

    if :blayers in which_plots
        xu = xs[idx_u[1]:idx_u[2]]
        xl = xs[idx_l[1]:-1:idx_l[2]]

        x_begin = 0.0
        if !isnothing(panels_wake)
            x_end = 2.0
        else
            x_end = 1.0
        end

        figure()
        ax1 = subplot(2,2,1)
        plot(xu, deltas_u, color=clr1, linewidth=2, label="upper")
        plot(xl, deltas_l, color=clr6, linewidth=2, label="lower")
        xlabel("normalized ξ")
        ylabel(L"\delta^*")
        xlim([x_begin, x_end])
        ax1.spines["right"].set_visible(false)
        ax1.spines["top"].set_visible(false)
        if isnothing(panels_wake); legend(); end

        ax2 = subplot(2,2,2)
        plot(xu, thetas_u, color=clr1, linewidth=2, label="upper")
        plot(xl, thetas_l, color=clr6, linewidth=2, label="lower")
        xlabel("normalized ξ")        
        ylabel(L"\theta^*")
        xlim([x_begin, x_end])
        ax2.spines["right"].set_visible(false)
        ax2.spines["top"].set_visible(false)

        ax3 = subplot(2,2,3)
        plot(xu[1:idx_transition_u-1], state3s_u[1:idx_transition_u-1], color=clr1, linewidth=2, label="upper")
        plot(xl[1:idx_transition_l-1], state3s_l[1:idx_transition_l-1], color=clr6, linewidth=2, label="lower")
        xlabel("normalized ξ")
        ylabel(L"\tilde{n}")
        xlim([x_begin, x_end])
        ax3.spines["right"].set_visible(false)
        ax3.spines["top"].set_visible(false)

        ax4 = subplot(2,2,4)
        if stu != 0.0; plot(xu[idx_transition_u:end], state3s_u[idx_transition_u:end] .^ 2, color=clr1, linewidth=2, label="upper"); end
        if stl != 0.0; plot(xl[idx_transition_l:end], state3s_l[idx_transition_l:end] .^ 2, color=clr6, linewidth=2, label="lower"); end
        xlabel("normalized ξ")
        ylabel(L"c_\tau")
        xlim([x_begin, x_end])
        ax4.spines["right"].set_visible(false)
        ax4.spines["top"].set_visible(false)

        if !isnothing(panels_wake)
            ax1.plot(xs_wake, deltas_w, color=clr5, linewidth=2, label="wake")
            ax1.legend()

            ax2.plot(xs_wake, thetas_w, color=clr5, linewidth=2)

            ax4.plot(xs_wake, state3s_w .^ 2, color=clr5, linewidth=2)
        end

        if data_file != nothing
            xu_mfoil = x_mfoil[idx_u_mfoil]
            xl_mfoil = x_mfoil[idx_l_mfoil]
            xw_mfoil = x_mfoil[idx_w_mfoil]

            ax1.plot(xu_mfoil, deltas_mfoil[idx_u_mfoil], color=clr1, "--")
            ax1.plot(xl_mfoil, deltas_mfoil[idx_l_mfoil], color=clr6, "--")
            ax1.plot(xw_mfoil, deltas_mfoil[idx_w_mfoil], color=clr5, "--")

            ax2.plot(xu_mfoil, thetas_mfoil[idx_u_mfoil], color=clr1, "--")
            ax2.plot(xl_mfoil, thetas_mfoil[idx_l_mfoil], color=clr6, "--")
            ax2.plot(xw_mfoil, thetas_mfoil[idx_w_mfoil], color=clr5, "--")

            x_3_u_mfoil = xu_mfoil[turb_mfoil[idx_u_mfoil] .== 0]
            ns_u_mfoil = state3s_mfoil[idx_u_mfoil][turb_mfoil[idx_u_mfoil] .== 0]
            x_3_l_mfoil= xl_mfoil[turb_mfoil[idx_l_mfoil] .== 0]
            ns_l_mfoil = state3s_mfoil[idx_l_mfoil][turb_mfoil[idx_l_mfoil] .== 0]

            x_4_u_mfoil = xu_mfoil[turb_mfoil[idx_u_mfoil] .== 1]
            cts_u_mfoil = state3s_mfoil[idx_u_mfoil][turb_mfoil[idx_u_mfoil] .== 1]
            x_4_l_mfoil= xl_mfoil[turb_mfoil[idx_l_mfoil] .== 1]
            cts_l_mfoil = state3s_mfoil[idx_l_mfoil][turb_mfoil[idx_l_mfoil] .== 1]

            cts_w_mfoil = state3s_mfoil[idx_w_mfoil]

            # subplot(2,2,3)
            ax3.plot(x_3_u_mfoil, ns_u_mfoil, color=clr1, "--")
            ax3.plot(x_3_l_mfoil, ns_l_mfoil, color=clr6, "--")

            # subplot(2,2,4)
            ax4.plot(x_4_u_mfoil, cts_u_mfoil .^ 2, color=clr1, "--")
            ax4.plot(x_4_l_mfoil, cts_l_mfoil .^ 2, color=clr6, "--")
            ax4.plot(xw_mfoil, cts_w_mfoil .^ 2, color=clr5, "--")
        end

        tight_layout()
    end

    if :blayers_paper in which_plots

        xu = xs[idx_u[1]:idx_u[2]]
        xl = xs[idx_l[1]:-1:idx_l[2]]
        
        if !(:blayers in which_plots)
            x_begin = 0.0
            if !isnothing(panels_wake)
                x_end = 2.0
            else
                x_end = 1.0
            end
        end

        figure(figsize=(8,6))

        if data_file != nothing

            xu_mfoil = x_mfoil[idx_u_mfoil]
            xl_mfoil = x_mfoil[idx_l_mfoil]
            xw_mfoil = x_mfoil[idx_w_mfoil]

            ax1 = subplot(2,2,1)
            plot(xu_mfoil, deltas_mfoil[idx_u_mfoil], linewidth=2, color=clr1, label="mfoil")
            plot(xl_mfoil, deltas_mfoil[idx_l_mfoil], linewidth=2, color=clr1)
            plot(xw_mfoil, deltas_mfoil[idx_w_mfoil], linewidth=2, color=clr1)
            legend()
            ax1.spines["right"].set_visible(false)
            ax1.spines["top"].set_visible(false)

            ax2 = subplot(2,2,2)
            plot(xu_mfoil, thetas_mfoil[idx_u_mfoil], linewidth=2, color=clr1)
            plot(xl_mfoil, thetas_mfoil[idx_l_mfoil], linewidth=2, color=clr1)
            plot(xw_mfoil, thetas_mfoil[idx_w_mfoil], linewidth=2, color=clr1)
            ax2.spines["right"].set_visible(false)
            ax2.spines["top"].set_visible(false)

            x_3_u_mfoil = xu_mfoil[turb_mfoil[idx_u_mfoil] .== 0]
            ns_u_mfoil = state3s_mfoil[idx_u_mfoil][turb_mfoil[idx_u_mfoil] .== 0]
            x_3_l_mfoil= xl_mfoil[turb_mfoil[idx_l_mfoil] .== 0]
            ns_l_mfoil = state3s_mfoil[idx_l_mfoil][turb_mfoil[idx_l_mfoil] .== 0]

            x_4_u_mfoil = xu_mfoil[turb_mfoil[idx_u_mfoil] .== 1]
            cts_u_mfoil = state3s_mfoil[idx_u_mfoil][turb_mfoil[idx_u_mfoil] .== 1]
            x_4_l_mfoil= xl_mfoil[turb_mfoil[idx_l_mfoil] .== 1]
            cts_l_mfoil = state3s_mfoil[idx_l_mfoil][turb_mfoil[idx_l_mfoil] .== 1]

            cts_w_mfoil = state3s_mfoil[idx_w_mfoil]

            ax3 = subplot(2,2,3)
            plot(x_3_u_mfoil, ns_u_mfoil, linewidth=2, color=clr1)
            plot(x_3_l_mfoil, ns_l_mfoil, linewidth=2, color=clr1)
            ax3.spines["right"].set_visible(false)
            ax3.spines["top"].set_visible(false)

            ax4 = subplot(2,2,4)
            plot(x_4_u_mfoil, cts_u_mfoil .^ 2, linewidth=2, color=clr1)
            plot(x_4_l_mfoil, cts_l_mfoil .^ 2, linewidth=2, color=clr1)
            plot(xw_mfoil, cts_w_mfoil .^ 2, linewidth=2, color=clr1)
            ax4.spines["right"].set_visible(false)
            ax4.spines["top"].set_visible(false)
        end

        ax1 = subplot(2,2,1)
        plot(xu, deltas_u, color=clr6, "--", linewidth=2, label="FLOWFoil")
        plot(xl, deltas_l, color=clr6, "--", linewidth=2)
        xlabel("normalized ξ")
        ylabel(L"\delta^*")
        xlim([x_begin, x_end])
        legend()

        ax2 = subplot(2,2,2)
        plot(xu, thetas_u, color=clr6, "--", linewidth=2)
        plot(xl, thetas_l, color=clr6, "--", linewidth=2)
        xlabel("normalized ξ")        
        ylabel(L"\theta^*")
        xlim([x_begin, x_end])

        ax3 = subplot(2,2,3)
        plot(xu[1:idx_transition_u-1], state3s_u[1:idx_transition_u-1], color=clr6, "--", linewidth=2)
        plot(xl[1:idx_transition_l-1], state3s_l[1:idx_transition_l-1], color=clr6, "--", linewidth=2)
        xlabel("normalized ξ")
        ylabel(L"\tilde{n}")
        xlim([x_begin, x_end])

        ax4 = subplot(2,2,4)
        if stu != 0.0; plot(xu[idx_transition_u:end], state3s_u[idx_transition_u:end] .^ 2, color=clr6, "--", linewidth=2); end
        if stl != 0.0; plot(xl[idx_transition_l:end], state3s_l[idx_transition_l:end] .^ 2, color=clr6, "--", linewidth=2); end
        xlabel("normalized ξ")
        ylabel(L"c_\tau")
        xlim([x_begin, x_end])

        if !isnothing(panels_wake)
            ax1.plot(xs_wake, deltas_w, color=clr6, "--", linewidth=2)
            ax2.plot(xs_wake, thetas_w, color=clr6, "--", linewidth=2)
            ax4.plot(xs_wake, state3s_w .^ 2, color=clr6, "--", linewidth=2)
        end

        tight_layout()
    end

    if :sweep_alpha in which_plots #inviscidoutputs is a vector of inviscidoutputs for every alpha
        cls = inviscidoutputs.cl
        cds = isnothing(viscousoutputs) ? zeros(length(cls)) : viscousoutputs.cd
        cms = inviscidoutputs.cmqc
        alphas = alphas * 180/pi

        figure(figsize=[12,4])
        subplot(1,3,1)
        plot(alphas, cls, color=clr1, linewidth=2, label=L"c_l")
        ylabel(L"c_l")
        xlabel(L"\alpha\ (rad)")

        subplot(1,3,2)
        plot(alphas, cds, color=clr1, linewidth=2, label=L"c_d")
        ylabel(L"c_d")
        xlabel(L"\alpha\ (rad)")

        subplot(1,3,3)
        plot(alphas, cms, color=clr1, linewidth=2, label=L"c_{m_{1/4}}")
        ylabel(L"c_{m_{1/4}}")
        xlabel(L"\alpha\ (rad)")

        tight_layout()

        #todo add data
    end

    if :sweep_Re in which_plots #inviscidoutputs is a vector of inviscidoutputs for every alpha
        cls = [o.cl for o in inviscidoutputs]
        cds = [o.cd for o in inviscidoutputs]
        cms = [o.cmqc for o in inviscidoutputs]

        figure()
        subplot(1,3,1)
        plot(Res, cls, color=clr1, linewidth=2, label=L"c_l")
        ylabel(L"c_l")
        xlabel(L"Re")

        subplot(1,3,2)
        plot(Res, cds, color=clr1, linewidth=2, label=L"c_d")
        ylabel(L"c_d")
        xlabel(L"Re")

        subplot(1,3,3)
        plot(Res, cms, color=clr1, linewidth=2, label=L"c_{m_{1/4}}")
        ylabel(L"c_{m_{1/4}}")
        xlabel(L"Re")
    end

    if :sweep_Mach in which_plots #inviscidoutputs is a vector of inviscidoutputs for every alpha
        cls = [o.cl for o in inviscidoutputs]
        cds = [o.cd for o in inviscidoutputs]
        cms = [o.cmqc for o in inviscidoutputs]

        figure()
        subplot(1,3,1)
        plot(Machs, cls, color=clr1, linewidth=2, label=L"c_l")
        ylabel(L"c_l")
        xlabel(L"Mach")

        subplot(1,3,2)
        plot(Machs, cds, color=clr1, linewidth=2, label=L"c_d")
        ylabel(L"c_d")
        xlabel(L"Mach")

        subplot(1,3,3)
        plot(Machs, cms, color=clr1, linewidth=2, label=L"c_{m_{1/4}}")
        ylabel(L"c_{m_{1/4}}")
        xlabel(L"Mach")
    end

    if :sweep_rho in which_plots #inviscidoutputs is a vector of inviscidoutputs for every alpha
        cls = [o.cl for o in inviscidoutputs]
        cds = [o.cd for o in inviscidoutputs]
        cms = [o.cmqc for o in inviscidoutputs]

        figure()
        subplot(1,3,1)
        plot(rhos, cls, color=clr1, linewidth=2, label=L"c_l")
        ylabel(L"c_l")
        xlabel(L"\rho")

        subplot(1,3,2)
        plot(rhos, cds, color=clr1, linewidth=2, label=L"c_d")
        ylabel(L"c_d")
        xlabel(L"\rho")

        subplot(1,3,3)
        plot(rhos, cms, color=clr1, linewidth=2, label=L"c_{m_{1/4}}")
        ylabel(L"c_{m_{1/4}}")
        xlabel(L"\rho")
    end

    if :sweep_Vinf in which_plots #inviscidoutputs is a vector of inviscidoutputs for every alpha
        cls = [o.cl for o in inviscidoutputs]
        cds = [o.cd for o in inviscidoutputs]
        cms = [o.cmqc for o in inviscidoutputs]

        figure()
        subplot(1,3,1)
        plot(Vinfs, cls, color=clr1, linewidth=2, label=L"c_l")
        ylabel(L"c_l")
        xlabel(L"V_\infty")

        subplot(1,3,2)
        plot(Vinfs, cds, color=clr1, linewidth=2, label=L"c_d")
        ylabel(L"c_d")
        xlabel(L"V_\infty")

        subplot(1,3,3)
        plot(Vinfs, cms, color=clr1, linewidth=2, label=L"c_{m_{1/4}}")
        ylabel(L"c_{m_{1/4}}")
        xlabel(L"V_\infty")
    end

    # if verbose
    #     println("cl: $(inviscidoutputs.cl)")
    #     println("cd: $(inviscidoutputs.cd)")
    #     println("cmqc: $(inviscidoutputs.cmqc)")
    # end
end

function sweep_run_inviscid_alpha(airfoilfile::String, numpanels::Int64, Re, Mach, rho, Vinf, chord, alphas::AbstractVector; name="")

    # get Airfoil and Panels objects
    airfoil = get_airfoil_from_file(airfoilfile; name)
    panels = get_panels(airfoil, numpanels)

    operatingparameters = OperatingParameters(airfoilfile, numpanels, airfoil.gapTE, Re, Mach, rho, Vinf, chord, alphas, sin.(alphas), cos.(alphas), XfoilTyler.XFOIL(), false)

    outputs = inviscid.(operatingparameters, Ref(panels))

    plotting_functions(panels, outputs, [:sweep_alpha]; data_file, verbose=false, alphas)
end

function sweep_run_inviscid_Re(airfoilfile::String, numpanels::Int64, Res::AbstractVector, Mach, rho, Vinf, chord, alpha; name="")

    operatingparameters = OperatingParameters.(Ref(airfoilfile), Ref(numpanels), Ref(false), Re, Ref(Mach), Ref(rho), Ref(Vinf), Ref(chord), alpha, sin(alpha), cos(alpha))

    # get Airfoil and Panels objects
    airfoil = get_airfoil_from_file(operatingparameters[1].airfoilfile; name)
    panels = get_panels(airfoil, operatingparameters[1].numpanels)

    outputs = inviscid.(operatingparameters, Ref(panels))

    plotting_functions(panels, outputs, [:sweep_Re]; data_file, verbose=false, Res)
end

function sweep_run_inviscid_Mach(airfoilfile::String, numpanels::Int64, Re, Machs::AbstractVector, rho, Vinf, chord, alpha; name="")

    operatingparameters = OperatingParameters.(Ref(airfoilfile), Ref(numpanels), Ref(false), Ref(Re), Machs, Ref(rho), Ref(Vinf), Ref(chord), alpha, sin(alpha), cos(alpha))

    # get Airfoil and Panels objects
    airfoil = get_airfoil_from_file(operatingparameters[1].airfoilfile; name)
    panels = get_panels(airfoil, operatingparameters[1].numpanels)

    outputs = inviscid.(operatingparameters, Ref(panels))

    plotting_functions(panels, outputs, [:sweep_Mach]; data_file, verbose=false, Machs)
end

function sweep_run_inviscid_rho(airfoilfile::String, numpanels::Int64, Re, Mach, rhos::AbstractVector, Vinf, chord, alpha; name="")

    operatingparameters = OperatingParameters.(Ref(airfoilfile), Ref(numpanels), Ref(false), Ref(Re), Ref(Mach), rhos, Ref(Vinf), Ref(chord), alpha, sin(alpha), cos(alpha))

    # get Airfoil and Panels objects
    airfoil = get_airfoil_from_file(operatingparameters[1].airfoilfile; name)
    panels = get_panels(airfoil, operatingparameters[1].numpanels)

    outputs = inviscid.(operatingparameters, Ref(panels))

    plotting_functions(panels, outputs, [:sweep_rho]; data_file, verbose=false, Machs)
end

function sweep_run_inviscid_Vinf(airfoilfile::String, numpanels::Int64, Re, Mach, rho, Vinfs::AbstractVector, chord, alpha; name="")

    operatingparameters = OperatingParameters.(Ref(airfoilfile), Ref(numpanels), Ref(false), Ref(Re), Ref(Mach), Ref(rho), Vinfs, Ref(chord), alpha, sin(alpha), cos(alpha))

    # get Airfoil and Panels objects
    airfoil = get_airfoil_from_file(operatingparameters[1].airfoilfile; name)
    panels = get_panels(airfoil, operatingparameters[1].numpanels)

    outputs = inviscid.(operatingparameters, Ref(panels))

    plotting_functions(panels, outputs, [:sweep_Vinf]; data_file, verbose=false, Machs)
end

