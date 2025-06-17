using Plots; pgfplotsx()
using LaTeXStrings

default()
default(
#     #:Plot
    # background_color = nothing,
    # background_color_outside = nothing,
#     display_type,
#     dpi,
#     extra_kwargs,
#     extra_plot_kwargs,
#     fontfamily,
#     foreground_color,
#     html_output_format,
#     inset_subplots,
#     layout,
#     link,
#     overwrite_figure,
#     plot_title,
#     plot_title_location,
#     plot_titlefontcolor,
#     plot_titlefontfamily,
#     plot_titlefonthalign,
#     plot_titlefontrotation,
#     plot_titlefontsize,
#     plot_titlefontvalign,
#     pos,
#     show,
    size = (200,200),
#     tex_output_standalone,
#     thickness_scaling,
#     warn_on_unsupported,
#     window_title,


#######################
#      :Series
#######################
#     arrow,
#     bar_edges,
#     bar_position,
#     bar_width,
#     bins,
#     colorbar_entry,
#     connections,
#     contour_labels,
#     contours,
#     extra_kwargs,
#     fill_z,
    fillalpha = 0.125,
    fillcolor = RGB(128/255,128/255,128/255),
#     fillrange,
#     group,
#     hover,
#     label,
#     levels,
#     line_z,
#     linealpha,
#     linecolor,
#     linestyle,
    linewidth = 2.0,
#     marker_z,
#     markeralpha,
#     markercolor,
#     markershape,
#     markersize,
    markerstrokealpha = 0,
#     markerstrokecolor,
#     markerstrokestyle,
#     markerstrokewidth,
#     normalize,
#     orientation,
#     primary,
#     quiver,
#     ribbon,
#     series_annotations,
#     seriesalpha,
#     seriescolor,
#     seriestype,
#     show_empty_bins,
#     smooth,
#     stride,
#     subplot,
#     weights,
#     x,
#     xerror,
#     y,
#     yerror,
#     y,
#     zerror


#######################
#      :Subplot
#######################
#     annotationcolor,
#     annotationfontfamily,
#     annotationfontsize,
#     annotationhalign,
#     annotationrotation,
#     annotations,
#     annotationvalign,
#     aspect_ratio,
    background_color_inside = nothing,
    background_color_legend = nothing,
    background_color_subplot = nothing,
#     bottom_margin,
#     camera,
#     clims,
    color_palette = [RGB(0,46/255,93/255), RGB(155/255,0,0)],
#     colorbar,
#     colorbar_continuous_values,
#     colorbar_discrete_values,
#     colorbar_fontfamily,
#     colorbar_formatter,
#     colorbar_scale,
#     colorbar_tickfontcolor,
#     colorbar_tickfontfamily,
#     colorbar_tickfonthalign,
#     colorbar_tickfontrotation,
#     colorbar_tickfontsize,
#     colorbar_tickfontvalign,
#     colorbar_ticks,
#     colorbar_title,
#     colorbar_title_location,
#     colorbar_titlefontcolor,
#     colorbar_titlefontfamily,
#     colorbar_titlefonthalign,
#     colorbar_titlefontrotation,
#     colorbar_titlefontsize,
#     colorbar_titlefontvalign,
#     extra_kwargs,
#     fontfamily_subplot,
    foreground_color_legend = nothing,
#     foreground_color_subplot,
#     foreground_color_title,
#     framestyle = :zerolines,
#     left_margin,
    legend = false, # include legend true/false
#     legendfontcolor,
#     legendfontfamily,
#     legendfonthalign,
#     legendfontrotation,
#     legendfontsize,
#     legendfontvalign,
#     legendtitle,
#     legendtitlefontcolor,
#     legendtitlefontfamily,
#     legendtitlefonthalign,
#     legendtitlefontrotation,
#     legendtitlefontsize,
#     legendtitlefontvalign,
#     margin,
#     projection,
#     right_margin,
#     subplot_index,
#     title,
#     titlefontcolor,
#     titlefontfamily,
#     titlefonthalign,
#     titlefontrotation,
#     titlefontsize,
#     titlefontvalign,
#     titlelocation,
#     top_margin


#####################
#       :Axis
#####################
#     discrete_values,
#     draw_arrow,
#     flip,
#     foreground_color_axis,
#     foreground_color_border,
#     foreground_color_grid,
#     foreground_color_guide,
#     foreground_color_minor_grid,
#     foreground_color_text,
#     formatter,
    grid = false, # background grid true/false
#     gridalpha,
    gridlinewidth = 0.5,
#     gridstyle,
#     guide,
#     guide_position,
#     guidefontcolor,
#     guidefontfamily,
#     guidefonthalign,
#     guidefontrotation,
#     guidefontsize,
#     guidefontvalign,
    ylims=(0,3),
    xlims=(0,2),
#     link,
#     minorgrid,
#     minorgridalpha,
#     minorgridlinewidth,
#     minorgridstyle,
#     minorticks,
#     mirror,
#     rotation,
    # scale,
    # showaxis = false, #turns off spines and tick labels, but not ticks
#     tick_direction,
#     tickfontcolor,
#     tickfontfamily,
#     tickfonthalign,
#     tickfontrotation,
#     tickfontsize,
#     tickfontvalign,
    ticks = false, #turns off tick marks
#     widen,
    )


## -- Constant Source Distribution Plot
x = [0;1.75]
y = [1.25;1.25]

plot1 = Plots.plot(x,y,fillrange=[0;0],fillcolor=1,fillalpha=0.125,xlabel="s",ylabel=L"$q(s)=$ const");

savefig(plot1, "potentialflowcontents/potentialflowfigures/constantsourcedist.tikz")


## -- Linear Source Distribution Plot

x = [0;1.75]
y = [1.25;1.25]

plot2 = Plots.plot(x,y,fillrange=[0;0],fillcolor=1,xlabel="s",ylabel=L"$q(s)=$ linear");

y2 = [1.25;2.25]
Plots.plot!(x,y2,fillrange=[1.25;1.25],fillcolor=2,);

savefig(plot2, "potentialflowcontents/potentialflowfigures/linearsourcedist.tikz")



## -- Uniform Flow

x = [0;2]
y = [0.05;0.25]


uplot = Plots.plot(x,y,linecolor=1,showaxis=false,arrow=true)
for i=1:7
    Plots.plot!(x,(i*0.35).+y,linecolor=1,arrow=true)
end

savefig(uplot, "potentialflowcontents/potentialflowfigures/uniformflow.tikz")


## -- Source Flow
function makecircle(R,h,v,n)
    t = range(pi/2,stop=5*pi/2,length=n)
    x = R*cos.(t) .+ h
    y = R*sin.(t) .+ v

    return x,y
end

x,y = makecircle(1.5,1.5,1.5,12)

s1plot = Plots.plot([1.5],[1.5],markershape=:circle,showaxis=false,aspect_ratio=:square,xlims=(0,3));
for i=1:12
    Plots.plot!([1.5;x[i]],[1.5,y[i]],linecolor=1,arrow=true,aspect_ratio=:square)
end

savefig(s1plot, "potentialflowcontents/potentialflowfigures/sourceflow.tikz")

## -- Sink Flow
x1,y1 = makecircle(1.5,1.5,1.5,12)
x2,y2 = makecircle(0.35,1.5,1.5,12)

s1plot = Plots.plot([1.5],[1.5],markershape=:circle,showaxis=false,aspect_ratio=:square,xlims=(0,3));
for i=1:12
    Plots.plot!([x1[i];x2[i]],[y1[i],y2[i]],linecolor=1,arrow=true,aspect_ratio=:square)
end

savefig(s1plot, "potentialflowcontents/potentialflowfigures/sinkflow.tikz")


## -- Vortex Flow

vplot = Plots.plot([1.5],[1.5],markershape=:circle,showaxis=false,aspect_ratio=:square,xlims=(0,3));
for i=1:4
    x,y = makecircle(i*0.35,1.5,1.5,120)
    Plots.plot!(x[1:end-2],y[1:end-2],linecolor=1,arrow=true,aspect_ratio=:square)
end

savefig(vplot, "potentialflowcontents/potentialflowfigures/vortexflow.tikz")



## -- doublet Flow

dplot = Plots.plot([0],[0],markershape=:circle,showaxis=false,aspect_ratio=:square,xlims=(-1.5,1.5),ylims=(-1.5,1.5));
for i=1:4
    x,y = makecircle(i*0.15,0,i*0.15,36)
    Plots.plot!(-x[1:end-1],y[1:end-1],linecolor=1,arrow=true,aspect_ratio=:square)
    Plots.plot!(-x[1:end-1],-y[1:end-1],linecolor=1,arrow=true,aspect_ratio=:square)
end

savefig(dplot, "potentialflowcontents/potentialflowfigures/doubletflow.tikz")