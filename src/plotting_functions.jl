############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Plotting gist geometry

function geo_plot(; 
        linewidth=2, 
        linewidth_axes=2, 
        saveplot=false,
        backgroundcolor = :grey95,
        config_name="config_name")

    xylabelsize=30
    xyticklabelsize=20

    f = Figure(size=(1600,1000))

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    global ax1 = Axis(f[1, 1], 
               title = latexstring(L"g^{xx} \textrm{ and } g^{xy}, \; \hat{s} = ", round(shat, digits=2)), 
               xlabel = L"θ / π",
               ylabel = L"a.u.",  
               titlesize = 30, 
               xlabelsize = xylabelsize, 
               ylabelsize = xylabelsize,
               xticklabelsize = xyticklabelsize,
               yticklabelsize = xyticklabelsize,
               backgroundcolor = backgroundcolor,
               xgridvisible = :false, 
               ygridvisible = :false) 

    p = lines!(th/pi, geo.gxx, linewidth = linewidth, color = :black, linestyle = :dash)
    p = lines!(th/pi, geo.gxy, linewidth = linewidth, color = :red, label = "gxy")
    p = hlines!(0, linewidth = linewidth_axes, color = :black)
    p = vlines!(0, linewidth = linewidth_axes, color = :black)

    hlines!(0, color = :black)

    global ax2 = Axis(f[2, 1], 
               title = latexstring(L"\textrm{Local shear: } \frac{d}{dθ} (g^{xy}/g^{xx})"), 
               xlabel = L"θ / π",
               ylabel = L"a.u.",  
               titlesize = 30, 
               xlabelsize = xylabelsize, 
               ylabelsize = xylabelsize,
               xticklabelsize = xyticklabelsize,
               yticklabelsize = xyticklabelsize,
               backgroundcolor = backgroundcolor,
               xgridvisible = :false, 
               ygridvisible = :false) 

    p = lines!(th/pi, geo.shat_loc, linewidth = linewidth, color = :blue)
    p = hlines!(0, linewidth = linewidth_axes, color = :black)
    p = vlines!(0, linewidth = linewidth_axes, color = :black)
    p = hlines!(shat, linewidth = linewidth, linestyle = :dash, color = :grey)

    hlines!(0, color = :black)

    global ax3 = Axis(f[1, 2], 
               title = latexstring(L"\textrm{∇B and curvature: }  \mathcal{L}_y \textrm{ and } \mathcal{K}_y"), 
               xlabel = L"θ / π",
               ylabel = L"a.u.",  
               titlesize = 30, 
               xlabelsize = xylabelsize, 
               ylabelsize = xylabelsize,
               xticklabelsize = xyticklabelsize,
               yticklabelsize = xyticklabelsize,
               backgroundcolor = backgroundcolor,
               xgridvisible = :false,
               ygridvisible = :false)

    p = lines!(th/pi, geo.Ly, linewidth = linewidth, color = :black, linestyle = :dash)
    p = lines!(th/pi, geo.Ky, linewidth = linewidth, color = :Green)
    p = hlines!(0, linewidth = linewidth_axes, color = :black)
    p = vlines!(0, linewidth = linewidth_axes, color = :black)

    global ax4 = Axis(f[2, 2], 
               title = latexstring(L"\textrm{Shear/FLR term: } g^{yy}"), 
               xlabel = L"θ / π",
               ylabel = L"a.u.",  
               titlesize = 30, 
               xlabelsize = xylabelsize, 
               ylabelsize = xylabelsize,
               xticklabelsize = xyticklabelsize,
               yticklabelsize = xyticklabelsize,
               backgroundcolor = backgroundcolor,
               xgridvisible = :false,
               ygridvisible = :false)

    p = lines!(th/pi, geo.gyy, linewidth = linewidth, color = :purple)
    p = hlines!(0, linewidth = linewidth_axes, color = :black)
    p = vlines!(0, linewidth = linewidth_axes, color = :black)

    global ax5 = Axis(f[1, 3], 
                      title = latexstring(L"\textrm{Magnetic field:} \hat{B}"), 
               xlabel = L"θ / π",
               ylabel = L"a.u.",  
               titlesize = 30, 
               xlabelsize = xylabelsize, 
               ylabelsize = xylabelsize,
               xticklabelsize = xyticklabelsize,
               yticklabelsize = xyticklabelsize,
               backgroundcolor = backgroundcolor,
               xgridvisible = :false,
               ygridvisible = :false)

    p = lines!(th/pi, geo.B, linewidth = linewidth, color = :turquoise2)
    p = hlines!(0, linewidth = linewidth_axes, color = :black)
    p = vlines!(0, linewidth = linewidth_axes, color = :black)

    global ax6 = Axis(f[2, 3], 
                      title = latexstring(L"\textrm{Jacobian:} \hat{J} = \sqrt{g}"), 
               xlabel = L"θ / π",
               ylabel = L"a.u.",  
               titlesize = 30, 
               xlabelsize = xylabelsize, 
               ylabelsize = xylabelsize,
               xticklabelsize = xyticklabelsize,
               yticklabelsize = xyticklabelsize,
               backgroundcolor = backgroundcolor,
               xgridvisible = :false,
               ygridvisible = :false)

    p = lines!(th/pi, geo.J, linewidth = linewidth, color = :Magenta)
    p = hlines!(0, linewidth = linewidth_axes, color = :black)
    p = vlines!(0, linewidth = linewidth_axes, color = :black)

    if saveplot==true
        save(string("./output/plots/geometry/",configname,"_global_shear_",round(shat, digits=2),"_dpdx_",round(my_dpdx, digits=2),"_geometry.png"),f)
    end

    f

end

# plotting array of eigenfunctions

function efn_array_plot(
        efn_matrix_input,
        x_var; 
        title=latexstring(L"\textrm{Title}"),
        xlabel=L"\textrm{x-label}",
        ylabel=L"\textrm{y-label}",
        text1=L"text1",
        text2=L"text2",
        text3=L"text3",
        text4=L"text4",
        text5=L"text5",
        text6=L"text6",
        text_x_position=1,
        text_y_position=1,
        col_select=[1,2,3,4,5,6],
        row_min_max=[1 2], 
        col_min_max=[1 3], 
        data_color=:red,
        textcolor=:red,
        textsize=30,
        linewidth=2, 
        linewidth_axes=2,
        linestyle=:solid,
        norm=true,
        saveplot=false,
        savefolder="./quick_plot.pdf")

    global textar = [text1, text2, text3, text4, text5, text6]

    f = Figure(size=(1600,1000))

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    for i in row_min_max[1]:row_min_max[2]
        
        for j in col_min_max[1]:col_min_max[2]

            global ax = Axis(f[i, j],
                             title=title,
                             xlabel=xlabel,
                             ylabel=ylabel,  
                             titlesize = 30, 
                             xlabelsize = 40, 
                             ylabelsize = 40,
                             xticklabelsize = 30,
                             yticklabelsize = 30,
                             xgridvisible = :false, 
                             ygridvisible = :false) 

            if norm==true
                p = lines!(x_var, efn_matrix_input[:,col_select[(i-1)*col_min_max[2] + j]] / maximum( efn_matrix_input[:,col_select[(i-1)*col_min_max[2] + j]]), linewidth=linewidth, linestyle=linestyle, color=data_color)
            elseif norm==false
                p = lines!(x_var, efn_matrix_input[:,col_select[(i-1)*col_min_max[2] + j]] / maximum( efn_matrix_input[:,col_select[(i-1)*col_min_max[2] + j]]), linewidth=linewidth, linestyle=linestyle, color=data_color)
            end

            t = text!(text_x_position, text_y_position, text = textar[(i-1)*col_min_max[2] + j], fontsize = textsize, color = textcolor)

            p = hlines!(0, linewidth = linewidth_axes, color = :black)
            p = vlines!(0, linewidth = linewidth_axes, color = :black)
    

            hlines!(0, color = :black)

        end

    end

    if saveplot==true
        save(savefolder,f; backend=CairoMakie)
    end

    f

end

function geo_plot_interactive(; 
        linewidth=3, 
        linewidth_axes=2, 
        linestyle=:solid, 
        linecolor=:blue)

    fig = Figure(size=(1000,700), font="Times New Roman")

    fontsize_theme = Theme(fontsize = 25)
    set_theme!(fontsize_theme)

    geo_data = [geo.gxx, geo.gxy, geo.gyy, geo.B, geo.J, geo.Ly, geo.Lx, geo.dBdz, geo.Ky, geo.shat_loc]
    geo_name = ["gxx", "gxy", "gyy", "B", "J", "Ly", "Lx", "dB/dz", "Ky", "Local shear"] 

    menu = Menu(fig,
                options = zip(geo_name, geo_data),
                default = "gyy")

    fig[1,1] = vgrid!(
                    Label(fig, L"\textrm{Geometric quantity}", width = nothing),
                    menu;
                    tellheight = false, width=220)

    global ax = Axis(fig[1, 2], 
               xlabel = L"θ / π",
               ylabel = L"a.u.",  
               titlesize = 30, 
               xlabelsize = 40, 
               ylabelsize = 40,
               xticklabelsize = 30,
               yticklabelsize = 30,
               backgroundcolor = :grey95,
               xgridvisible = :false, 
               ygridvisible = :false)

    data = Observable{Any}(geo_data[1])

    y = lift(data) do f
        f
    end

    lines!(th/pi, y, linewidth=linewidth, linestyle=linestyle, color=linecolor)
    hlines!(0, linewidth=linewidth_axes, linestyle=linestyle, color=:black)
    vlines!(0, linewidth=linewidth_axes, linestyle=linestyle, color=:black)

    on(menu.selection) do s
        data[] = s
        autolimits!(ax)
    end
    notify(menu.selection)

    fig

end

##################################### End ######################################
