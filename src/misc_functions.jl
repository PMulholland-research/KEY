############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Defining heaviside function
# H(x) = 1 if x >/= 0,
# H(x) = 0 if x < 0,

function heaviside_to_array(in_array, out_array)
    for i in 1:length(in_array)
        if in_array[i] >= 0
            append!(out_array, 1)
        else
            append!(out_array, 0)
        end
    end
end

function heaviside(in_array)
    for i in 1:length(in_array)
        if in_array[i] >= 0
            return(1)
        else
            return(0)
        end
    end
end

# Normalising input array to maximum value

function norm_max(input_ar)
    return input_ar ./ maximum(input_ar)
end

# Integrates array and appends result to an output array 
function integrate(interval, input_ar, output_ar; low_lim=1, up_lim=length(interval))
    append!(output_ar, trapz(interval[low_lim:up_lim], input_ar[low_lim:up_lim]))
end

function integrate_to_array(interval, input_ar, output_ar)
    for i in 1:length(interval)
        integrate(interval, input_ar, output_ar, up_lim=i)
    end
end

# Same as "integrate" function above, but doesn't require output array, just returns result
function qck_integrate(input_ar, interval; low_lim=1, up_lim=length(interval))
    return(trapz(interval[low_lim:up_lim], input_ar[low_lim:up_lim]))
end

function local_max(input_matrix)

    input_matrix_row_down = circshift(input_matrix, (1,0))
    input_matrix_row_up = circshift(input_matrix, (-1,0))
    input_matrix_col_right = circshift(input_matrix, (0,1))
    input_matrix_col_left = circshift(input_matrix, (0,-1))

    global local_max_index_mat = zeros(size(input_matrix)[1],size(input_matrix)[2])

    for i in 1:size(input_matrix)[1], j in 1:size(input_matrix)[2]

        if input_matrix[i,j] >= input_matrix_row_down[i,j] && input_matrix[i,j] >= input_matrix_row_up[i,j] && input_matrix[i,j] >= input_matrix_col_right[i,j] && input_matrix[i,j] >= input_matrix_col_left[i,j]
            local_max_index_mat[i,j] +=1
            #print("\nMatrix has local maximum at index (",i,",",j,") with value ",input_matrix[i,j],"\n")
        end
    end

    return local_max_index_mat
end

# Defining functions for Gauss-Legendre/Hermite/Laguerre

# Example: f(x) = x^4 # defining integrand function without exp term (exp is provided by quadrature scheme)

# Example: Gauss-Legendre quadrature
# ∫dx f(x) with limits [-1,1]
# check integral: with f(x) = x^4
# I = 2/5 = 0.4

function gauss_legendre(n)

    x, w = gausslegendre(n)

    I = dot(w, f.(x))

    print("\nGauss-Legendre integral result = ", I, "\n")

end

# Example: Gauss-Laguerre quadrature
# ∫dx exp(-x) f(x) with limits [0,+∞]
# check integral: with f(x) = x^4
# I = 24 

function gauss_laguerre(n)

    x, w = gausslaguerre(n)

    I = dot(w, f.(x))

    print("\nGauss-Laguerre integral result = ", I, "\n")

end

# Example: Gauss-Hermite quadrature
# ∫dx exp(-x^2) f(x) with limits [-∞,+∞]
# check integral: with f(x) = x^4
# I = 3√π/4 ≈ 1.329

function gauss_hermite(n)

    x, w = gausshermite(n)

    I = dot(w, f.(x))

    print("\nGauss-Hermite integral result = ", I, "\n")

end

function read_complex_datafile(filename)
    readdlm(filename, '\t', Complex{Float64})
end

# For reading eigenvalues.dat

function read_csv_dataframe(filename)
    CSV.read(filename, DataFrame) 
end

# For reading eigenvectors.dat

function read_complex_csv_dataframe(filename)
    CSV.read(filename, DataFrame, types=Complex{Float64}) 
end

# Quick-plot function (for output scanlogs, eigenvectors, geometry, etc.)

function quick_plot(x_data, y_data;
        x_label=L"\textrm{x-label}",
        y_label=L"\textrm{y-label}",
        x_tick_step=round(0.1*maximum(x_data), digits=1),
        y_tick_step=round(0.1*maximum(y_data), digits=1),
        x_ticks=0:x_tick_step:1.5*maximum(x_data),
        y_ticks=0:y_tick_step:1.5*maximum(y_data),
        add_text=false,
        x_text_pos=0.8*maximum(x_data),
        y_text_pos=0.8*maximum(y_data),
        plot_text=L"\textrm{Insert text here}",
        data_color=:red,
        linewidth=3,
        markersize=30,
        fontsize=50
    )

    global fig = Figure(size=(1200,800))

    ax = Axis(fig[1,1],
              xlabel=x_label,
              ylabel=y_label,
              xlabelsize=50,
              ylabelsize=50,
              xticks=x_ticks,
              yticks=y_ticks,
              xticklabelsize=30,
              yticklabelsize=30,
              xgridvisible=false,
              ygridvisible=false)

    lines!(ax, 
           x_data, 
           y_data, 
           color=data_color, 
           linewidth=linewidth)

    scatter!(ax, 
             x_data, 
             y_data, 
             color=data_color, 
             markersize=markersize)

    if add_text==true

        text!(ax, 
              x_text_pos, 
              y_text_pos, 
              text=plot_text, 
              color=data_color, 
              fontsize=fontsize)

    end

    hlines!(ax, 0, color=:black)
    vlines!(ax, 0, color=:black)

    xlims!(0,1.05*maximum(x_data))
    ylims!(0,1.05*maximum(y_data))

    fig

end

##################################### End ######################################
