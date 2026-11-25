############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Initialising nth-derivative arrays
function reset_ar()
    global proc_input_ar = Float64[] # initialising processed output array
    global proc_input_ar_extend = Float64[] # initialising processed output array
    global d1_ar = Float64[] # initialising 1st-derivative array
    global d2_ar = Float64[] # initialising 2nd-derivative array
    global xtra_d1_ar = Float64[] # initialising 1st-derivative array
    global xtra_d2_ar = Float64[] # initialising 2nd-derivative array
    global d1_mat_ar = Float64[] # initialising 1st-derivative array
    global d2_mat_ar = Float64[] # initialising 2nd-derivative array
end

reset_ar()

# Central finite difference (CFD) 
# Using stencils of various accuracy 
# See Wikipedia page: https://en.wikipedia.org/wiki/Finite_difference_coefficient

# Smoothing input data and appending additional points to start and end for CFD

function process_input_array(init_ar; 
        sf=0.1, # smoothing-intensity factor (needs to be Float64)
        norm_smooth=true, 
        just_smooth=false, 
        just_norm=false, 
        none=false,
        auto_arg=true, 
        extend=false,
        default_arg=th, 
        output_ar=proc_input_ar,
        output_ar_extend=proc_input_ar_extend)

    if norm_smooth == true
        global proc_ar = denoise(init_ar, factor=sf)[1] # smoothing input array (removes noise)
        global proc_ar = proc_ar ./ maximum(proc_ar) # normalising to 1
    elseif just_norm == true
        global proc_ar = init_ar ./ maximum(init_ar) # normalising to 1; no smoothing
    elseif just_smooth == true
        global proc_ar = denoise(init_ar, factor=sf)[1] # just smoothing; no normalising
    elseif none == true
        global proc_ar = 1*init_ar # no normalising, no smoothing (1 prefactor ensures input array isn't altered by processing!)
    end

    if auto_arg == true
        global init_arg = range(-pi,pi,length(init_ar))
    else
        global init_arg = default_arg
    end

    global hx = step(init_arg)

    append!(output_ar, proc_ar)
    global proc_arg_no_extend = range(first(init_arg),last(init_arg),length(proc_ar))

    if extend==true
        for i in 1:5
            pushfirst!(proc_ar, first(proc_ar))
            push!(proc_ar, last(proc_ar))
        end

        global proc_arg = range(first(init_arg)-(5*hx),last(init_arg)+(5*hx),length(proc_ar))

        append!(output_ar_extend, proc_ar)
    end

end

# First derivative

# CFD coefficient matrix for first derivatve
# row number x2: order of accuracy, e.g., 2nd row -> 4th order accuracy
# column number -6: nth term being considered, e.g., 8th column -> n=2 term, x_2
# x_n = x_0 + n*h_x; where n = (i-6) for these coefficient matrices below

CFD_coef_d1=[0 0 0 0 -1/2 0 1/2 0 0 0 0; 
             0 0 0 1/12 -2/3 0 2/3 -1/12 0 0 0;
             0 0 -1/60 3/20 -3/4 0 3/4 -3/20 1/60 0 0;
             0 1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280 0]

function loc_dfdx(; 
        proc_input_arg = proc_arg, 
        proc_input_ar = proc_ar, 
        x0=0, 
        acc=2, 
        output_ar = d1_ar) # x0: point where derivative is evaluated; acc: order of accuracy 

    f = LinearInterpolation(proc_input_arg, proc_input_ar)
    global loc = 0 # initialising 1st derivative at local point
    map(i -> loc += (1/hx) * CFD_coef_d1[Int64(acc/2),i]*f(x0 + (i-6)*hx), 1:11)
    push!(output_ar, loc)

end

function global_dfdx(; 
        init_input_arg = init_arg, 
        proc_input_arg = proc_arg, 
        proc_input_ar = proc_ar, 
        acc=2, 
        output_ar = d1_ar) 

    for i in 1:length(init_input_arg)
        loc_dfdx(x0 = init_input_arg[i], proc_input_arg=proc_input_arg, proc_input_ar=proc_input_ar, acc=acc, output_ar=output_ar)
    end

end

# Second derivative

CFD_coef_d2 =[0 0 0 0 1 -2 1 0 0 0 0; 
              0 0 0 -1/12 4/3 -5/2 4/3 -1/12 0 0 0;
              0 0 1/90 -3/20 3/2 -49/18 3/2 -3/20 1/90 0 0;
              0 -1/560 8/315 -1/5 8/5 -205/72 8/5 -1/5 8/315 -1/560 0]

function loc_d2fdx2(; 
        proc_input_arg = proc_arg, 
        proc_input_ar = proc_ar, 
        x0=0, # x0: point where derivative is evaluated 
        acc=2, # acc: order of accuracy
        output_ar = d2_ar) 

    f = LinearInterpolation(proc_input_arg, proc_input_ar)
    global loc = 0 # initialising 2nd derivative at local point
    map(i -> loc += (1/hx)^2 * CFD_coef_d2[Int64(acc/2),i]*f(x0 + (i-6)*hx), 1:11)
    push!(output_ar, loc)

end

function global_d2fdx2(; 
        init_input_arg = init_arg, 
        proc_input_arg = proc_arg, 
        proc_input_ar = proc_ar, 
        acc=2, 
        output_ar = d2_ar) 

    for i in 1:length(init_input_arg)
        loc_d2fdx2(x0 = init_input_arg[i], proc_input_arg=proc_input_arg, proc_input_ar=proc_input_ar, acc=acc, output_ar=output_ar)
    end

end

# Full-field line (global) derivatives using matrix multiplication 

# First derivative

function global_dfdx_matrix(;input_ar=zpar, output_ar=d1_mat_ar, acc=2)

    global d1_mat = zeros(length(th), length(th)+10); # initialising empty matrix for CFD 1st-derivative stencil
    map(i-> d1_mat[i,i.+(0:10)] = (1/hx) * CFD_coef_d1[Int64(acc/2),:], 1:length(th)); # filling CFD matrix with coefficients
    append!(output_ar, d1_mat*input_ar)

end

# Second derivative

function global_d2fdx2_matrix(;input_ar=zpar, output_ar=d2_mat_ar, acc=2)

    global d2_mat = zeros(length(th), length(th)+10); # initialising empty matrix for CFD 2nd-derivative stencil
    map(i-> d2_mat[i,i.+(0:10)] = (1/hx)^2 * CFD_coef_d2[Int64(acc/2),:], 1:length(th)); # filling CFD matrix with coefficients
    append!(output_ar, d2_mat*input_ar)

end

# Defining derivative matrices

function deriv_mat_cfd(;acc=2)

    global dfdx_matrix = zeros(length(th), length(th)+10); # initialising empty matrix for CFD 1st-derivative stencil
    map(i-> dfdx_matrix[i,i.+(0:10)] = (1/hx) * CFD_coef_d1[Int64(acc/2),:], 1:length(th));

    global d2fdx2_matrix = zeros(length(th), length(th)+10); # initialising empty matrix for CFD 2nd-derivative stencil
    map(i-> d2fdx2_matrix[i,i.+(0:10)] = (1/hx)^2 * CFD_coef_d2[Int64(acc/2),:], 1:length(th)); # filling CFD matrix with coefficients

end

# Forward/backward-finite-difference (F/BFD) stencils

# FFD

FFD_coef_d1=[-1 1 0 0 0 0 0 0 0;
             -3/2 2 -1/2 0 0 0 0 0 0;
             -11/6 3 -3/2 1/3 0 0 0 0 0;
             -25/12 4 -3 4/3 -1/4 0 0 0 0;
             -137/60 5 -5 10/3 -5/4 1/5 0 0 0;
             -49/20 6 -15/2 20/3 -15/4 6/5 -1/6 0 0]

FFD_coef_d2=[1 -2 1 0 0 0 0 0 0;
             2 -5 4 -1 0 0 0 0 0;
             35/12 -26/3 19/2 -14/3 11/12 0 0 0 0;
             15/4 -77/6 107/6 -13 61/12 -5/6 0 0 0;
             203/45 -87/5 117/4 -254/9 33/2 -27/5 137/180 0 0;
             469/90 -223/10 879/20 -949/18 41 -201/10 1019/180 -7/10 0]

BFD_coef_d1=[0 0 0 0 0 0 0 -1 1;
             0 0 0 0 0 0 1/2 -2 3/2;
             0 0 0 0 0 -1/3 3/2 -3 11/6;
             0 0 0 0 1/4 -4/3 3 -4 25/12;
             0 0 0 -1/5 5/4 -10/3 5 -5 137/60;
             0 0 1/6 -6/5 15/4 -20/3 15/2 -6 49/20]

BFD_coef_d2=[0 0 0 0 0 0 1 -2 1;
             0 0 0 0 0 -1 4 -5 2;
             0 0 0 0 11/12 -14/3 19/2 -26/3 35/12;
             0 0 0 -5/6 61/12 -13 107/6 -77/6 15/4;
             0 0 137/180 -27/5 33/2 -254/9 117/4 -87/5 203/45;
             0 -7/10 1019/180 -201/10 41 -949/18 879/20 -223/10 469/90]

function deriv_mat_FCB(;acc=2) # derivative matrix combining foward-central-backward (FCB) finite difference stencils

    hx = step(th)

    # Only using FFD/BFD at first/last row, and CFD everywhere else

    FFD_acc=2
    CFD_acc=2
    BFD_acc=2

    global dfdx_matrix_FCB = zeros(length(th), length(th)); # initialising empty matrix for 1st-derivative stencil

    map(i-> dfdx_matrix_FCB[i,i.+(0:8)] = (1/hx) * FFD_coef_d1[Int64(FFD_acc),:], 1:1); # FFD for left boundary
    map(i-> dfdx_matrix_FCB[i,i.+(-1:1)] = (1/hx) * CFD_coef_d1[Int64(CFD_acc/2),5:7], 2:length(th)-1); # CFD for middle bulk
    map(i-> dfdx_matrix_FCB[i,i.+(-8:0)] = (1/hx) * BFD_coef_d1[Int64(BFD_acc),:], length(th)); # BFD for right boundary

    global d2fdx2_matrix_FCB = zeros(length(th), length(th)); # initialising empty matrix for 2nd-derivative stencil

    map(i-> d2fdx2_matrix_FCB[i,i.+(0:8)] = (1/hx)^2 * FFD_coef_d2[Int64(FFD_acc),:], 1:1); # FFD for left boundary
    map(i-> d2fdx2_matrix_FCB[i,i.+(-1:1)] = (1/hx)^2 * CFD_coef_d2[Int64(CFD_acc/2),5:7], 2:length(th)-1); # CFD for middle bulk
    map(i-> d2fdx2_matrix_FCB[i,i.+(-8:0)] = (1/hx)^2 * BFD_coef_d2[Int64(BFD_acc),:], length(th)); # BFD for right boundary

end

##################################### End ######################################
