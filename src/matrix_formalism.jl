############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Matrix operator M ≡ (L-R) for KBM eigenvalue problem: MΦ=0
M_op_fluid(Ω,τ,kyrho,η_i,η_e,omn,β) = Tridiagonal(L_op - Diagonal(RHS_fluid(Ω,τ,kyrho,η_i,η_e,omn,β))) 
M_op_analytic_resonant(Ω,τ,kyrho,η_i,η_e,omn,β) = Tridiagonal(L_op - Diagonal(RHS_analytic_resonant(Ω,τ,kyrho,η_i,η_e,omn,β))) 
M_op_numeric_resonant(Ω,τ,kyrho,η_i,η_e,omn,β) = Tridiagonal(L_op - Diagonal(RHS_numeric_resonant(Ω,τ,kyrho,η_i,η_e,omn,β))) 

# NOTE: boundary conditions for MΦ=0 are imposed below

function reset_eval()
    global eval_condition_M_op = ComplexF64[]
end

# Using condition number of matrix to find root of Det[M(Ω)]→0 (or δDet[M(Ω)]≈0). 
# For a singular matrix (Det[M] → 0), the condition number goes to infinity.
# The following function finds Ω that maximises cond[M(Ω)] (i.e., Ω that brings it as close to infinity as possible).

function condition_M_op(β;
        kyrho_gene=0.2, 
        omti=3.5, 
        omte=0, 
        omn=1, 
        τ=1, 
        gamma_lims=[0,5],
        omega_lims=[0,5],
        gamma_step=0.1, 
        omega_step=0.1, 
        linewidth=2,
        linewidth_crosshairs=3,
        markersize_crosshairs=50,
        plot=true,
        axislines=false,
        guideline=false,
        crosshairs=true,
        analytic_resonant=false,
        numeric_resonant=false,
        n_ev=1,
        gene_units=true,
        quiet_time=true,
        quiet_output=false,
        progress_bar=true)
    
    kyrho = kyrho_key(kyrho_gene)

    global beta_condition_M_op = β
    global tau_condition_M_op = τ
    global kyrho_condition_M_op = kyrho
    global omti_condition_M_op = omti
    global omte_condition_M_op = omte
    global omn_condition_M_op = omn
    global analytic_resonant_condition_M_op = analytic_resonant
    global numeric_resonant_condition_M_op = numeric_resonant
    global gene_units_condition_M_op = gene_units

    # gradient ratios
    η_i = omti/omn
    η_e = omte/omn

    # setting limits and stepsize of evalmat
    gmin, gmx = gamma_lims 
    omin, omx = omega_lims 

    gstep = gamma_step
    ostep = omega_step

    reset_eval()

    global gene_xo = omin:ostep:omx
    global gene_xg = gmin:gstep:gmx

    xo = -Ω.(gene_xo,kyrho,omn,τ)
    xg = -Ω.(gene_xg,kyrho,omn,τ)

    # choosing calculation array: with/without progress bar
    progress_bar_array=ProgressBar(1:length(xo)) 
    for_loop_array=1:length(xo)

    if progress_bar==true
        calc_array=progress_bar_array 
    else
        calc_array=for_loop_array 
    end

    # filling evalmat with chosen eigenvalues
    global evalmat = ComplexF64.(xo) .+ transpose(ComplexF64.(im*xg))

    global condition_mat = Complex.(zeros(length(xo), length(xg)))

    if analytic_resonant==true && numeric_resonant==true
        throw(ArgumentError("You have selected BOTH analytic-resonant and numeric-resonant options. Please select only ONE option, or set both to false for non-resonant option."))
    end

    if analytic_resonant==false && numeric_resonant==false 

        t1 = time()
        for i in calc_array, j in 1:length(xg)

            global M_op = M_op_fluid(evalmat[i,j], τ, kyrho, η_i, η_e, omn, β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating condition number of M_op
            condition_mat[i,j] = cond(M_op,Inf) 

        end
        t2 = time()

    elseif analytic_resonant==true

        t1 = time()
        for i in calc_array, j in 1:length(xg)

            global M_op = M_op_analytic_resonant(evalmat[i,j], τ, kyrho, η_i, η_e, omn, β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating condition number of M_op
            condition_mat[i,j] = cond(M_op,Inf) 

        end
        t2 = time()

    elseif numeric_resonant==true

        t1 = time()
        for i in calc_array, j in 1:length(xg)

            global M_op = M_op_numeric_resonant(evalmat[i,j], τ, kyrho, η_i, η_e, omn, β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating condition number of M_op
            condition_mat[i,j] = cond(M_op,Inf) 

        end
        t2 = time()

    end

    condition_mat_no_NaN = 1*condition_mat
    condition_mat_no_NaN[isnan.(condition_mat_no_NaN)] .= 0

    global abs_condition_mat = abs.(condition_mat_no_NaN)

    if gene_units==true
    
        map(i-> append!(eval_condition_M_op, reverse(evalmat, dims=2)[local_max(reverse(abs_condition_mat, dims=2)) .== 1][i]), 1:n_ev)

        global eval_condition_M_op_gene_units = round.(-ω_norm.(eval_condition_M_op,kyrho,omn,τ), digits=3) # converts KEY eigenvalues to GENE units

        if quiet_output==false
            print(string("\nkyrho_gene = ",kyrho_gene))
            print(string("\n\nβ = ",round(β*100,digits=3),"%"))
            print(string("\n\nω_*pi/2 = ",round(-0.5*ω_norm(Ω_p_i(τ,η_i),kyrho,omn,τ),digits=4)))
            print(string("\n\nEigenvalues: \n"))
            for i in 1:n_ev
                print(i,") ",eval_condition_M_op_gene_units[i], "\n")         
                print("-----------------\n")         
            end
        end

        if quiet_time==false
            print("\nTime elapsed for calculating condition number: ", round(t2-t1, digits=3), " seconds.\n\n\n")
        end

        if plot==true

            f = Figure(size=(1000,1000))

            fontsize_theme = Theme(fontsize = 25)
            set_theme!(fontsize_theme)

            ax1 = Axis(f[1,1], 
                       title = latexstring(L"\textrm{δCondition[M(Ω)]≈0}", L"\; - \; β = ", round(β*100,digits=3), " %"),
                       xlabel = L"γ/(c_\mathrm{s}/L_\mathrm{ref})",
                       ylabel = L"ω_\mathrm{r}/(c_\mathrm{s}/L_\mathrm{ref})",
                       titlesize = 30, 
                       xlabelsize = 40, 
                       ylabelsize = 40,
                       xticklabelsize = 30,
                       yticklabelsize = 30)

            global condition_contour=contourf!(gene_xg, gene_xo, log10.(transpose(abs_condition_mat)), colormap = :plasma, levels=15)

            if axislines==true
                hlines!(0, color = :white)
                vlines!(0, color = :white)
            end

            if guideline==true
                hlines!(-0.5*ω_norm(Ω_p_i(τ,η_i),kyrho,omn,τ), color=:green, linestyle=:dash, linewidth=linewidth_crosshairs) # shows ω_r of strongly-driven KBM
            end

            if crosshairs==true
                scatter!(imag(-ω_norm.(eval_condition_M_op,kyrho,omn,τ)), real(-ω_norm.(eval_condition_M_op,kyrho,omn,τ)), color=:green, marker='◇', markersize=markersize_crosshairs)
                #map(i-> hlines!(real(-ω_norm.(eval_condition_M_op,kyrho,omn,τ)), color = :red, linestyle=:dash, linewidth=linewidth_crosshairs), 1:n_ev)
                #map(i-> vlines!(imag(-ω_norm.(eval_condition_M_op,kyrho,omn,τ)), color = :red, linestyle=:dash, linewidth=linewidth_crosshairs), 1:n_ev)
            end

            Colorbar(f[1,2], condition_contour , width=25, label=L"10^n") # log_10 scale

            f

        end

    elseif gene_units==false

        map(i-> append!(eval_condition_M_op, reverse(evalmat, dims=2)[local_max(reverse(abs_condition_mat, dims=2)) .== 1][i]), 1:n_ev)
    
        print("\nEigenvalues: \n")
        for i in 1:n_ev
            print(i,") ",round(eval_condition_M_op[i], digits=3), "\n")
        end

        if quiet_time==false
            print("\nTime elapsed for calculating condition number:\n", round(t2-t1, digits=3), " seconds\n")
        end

        if plot==true

            f = Figure(size=(1000,1000))

            fontsize_theme = Theme(fontsize = 25)
            set_theme!(fontsize_theme)

            ax1 = Axis(f[1,1], 
                       title = latexstring(L"\textrm{δCondition[M(Ω)]≈0}", L"\; - \; β = ", round(β*100,digits=3), " %"),
                       xlabel = L"Ω_\mathrm{i} = γ/ω_\mathrm{*e}",
                       ylabel = L"Ω_\mathrm{r} = ω_r/ω_\mathrm{*e}",  
                       titlesize = 30, 
                       xlabelsize = 40, 
                       ylabelsize = 40,
                       xticklabelsize = 30,
                       yticklabelsize = 30)

            global condition_contour=contourf!(xg, xo, log10.(transpose(abs_condition_mat)), colormap = :plasma, levels=15)

            if axislines==true
                hlines!(0, color = :white)
                vlines!(0, color = :white)
            end

            if guideline==true
                hlines!(0.5*Ω_p_i(τ,η_i), color=:green, linestyle=:dash, linewidth=linewidth_crosshairs) # shows Ω_r of strongly-driven KBM
            end

            if crosshairs==true
                scatter!(imag(eval_condition_M_op), real(eval_condition_M_op), color=:green, marker='◇', markersize=markersize_crosshairs)
                #map(i-> hlines!(real(eval_condition_M_op[i]), color = :red, linestyle=:dash, linewidth=linewidth_crosshairs), 1:n_ev)
                #map(i-> vlines!(imag(eval_condition_M_op[i]), color = :red, linestyle=:dash, linewidth=linewidth_crosshairs), 1:n_ev)
            end

            Colorbar(f[1,2], condition_contour , width=25, label=L"10^n") # log_10 scale

            f

        end

    end

end

function reset_eval_evec_solutions_M_op()
    global eval_solutions_M_op = []
    global evec_solutions_M_op = []
end

# To be used after running condition_M_op() above 

function eval_evec_soln_condition_M_op(;
        β=beta_condition_M_op, 
        τ=tau_condition_M_op,
        kyrho=kyrho_condition_M_op,
        omti=omti_condition_M_op, 
        omte=omte_condition_M_op,
        omn=omn_condition_M_op,
        analytic_resonant=analytic_resonant_condition_M_op,
        numeric_resonant=numeric_resonant_condition_M_op,
        compute_all_evecs=true, # computing eigenvectors for all computed eigenvalues
        single_soln_index=1, # computing eigenvector for single specific eigenvalue solutiun
        min_soln=1, 
        max_soln=length(eval_condition_M_op),
        quiet_time=true)

    reset_eval_evec_solutions_M_op()

    η_i = omti/omn
    η_e = omte/omn
    
    if compute_all_evecs==true
        soln_range = min_soln:max_soln
    else
        soln_range = single_soln_index 
    end

    if analytic_resonant==false && numeric_resonant==false

        t1 = time()
        for i in soln_range 

            # Defining shorthand matrix: M_op
            global M_op = M_op_fluid(eval_condition_M_op[i],τ,kyrho,η_i,η_e,omn,β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating λ-eigenvalues and Φ eigenvectors
            push!(eval_solutions_M_op, eigvals(M_op))
            push!(evec_solutions_M_op, eigvecs(M_op))

        end
        t2 = time()

    elseif analytic_resonant==true && numeric_resonant==false

        t1 = time()
        for i in soln_range 

            # Defining shorthand matrix: M_op
            global M_op = M_op_analytic_resonant(eval_condition_M_op[i],τ,kyrho,η_i,η_e,omn,β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating λ-eigenvalues and Φ eigenvectors
            push!(eval_solutions_M_op, eigvals(M_op))
            push!(evec_solutions_M_op, eigvecs(M_op))

        end
        t2 = time()

    elseif numeric_resonant==true && analytic_resonant==false

        t1 = time()
        for i in soln_range 

            # Defining shorthand matrix: M_op
            global M_op = M_op_numeric_resonant(eval_condition_M_op[i],τ,kyrho,η_i,η_e,omn,β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating λ-eigenvalues and Φ eigenvectors
            push!(eval_solutions_M_op, eigvals(M_op))
            push!(evec_solutions_M_op, eigvecs(M_op))

        end
        t2 = time()

    end

    if quiet_time==false
        print("\nTime elapsed for calculating eigenvalues and eigenvectors:\n", round(t2-t1, digits=2), " seconds\n")
    end

end

# Can be run at any time: user inputs eigenvalue of interest to find corresponding eigenvector

function eval_evec_soln_M_op_any_input(eigenvalue; # must be in GENE units; accepts single eigenvalue and arrays of eigenvalues
        β=0.01, 
        kyrho_gene=0.2,
        omti=3.5, 
        omte=0,
        omn=1,
        τ=1,
        analytic_resonant=0,
        numeric_resonant=0,
        quiet_time=true)

    reset_eval_evec_solutions_M_op()

    kyrho = kyrho_key(kyrho_gene)

    η_i = omti/omn
    η_e = omte/omn
    
    # checking if input is single eigenvalue or array of eigenvalues
    if typeof(eigenvalue)==ComplexF64 || typeof(eigenvalue)==Complex{Int64}
        global eigenvalue_array = -Ω.(ComplexF64[eigenvalue],kyrho,omn,τ)
    elseif typeof(eigenvalue)==Vector{ComplexF64}
        global eigenvalue_array = -Ω.(eigenvalue,kyrho,omn,τ)
    end

    if analytic_resonant==false && numeric_resonant==false

        t1 = time()
        for i in 1:length(eigenvalue_array)

            # Defining shorthand matrix: M_op
            global M_op = M_op_fluid(eigenvalue_array[i],τ,kyrho,η_i,η_e,omn,β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating λ-eigenvalues and Φ eigenvectors
            push!(eval_solutions_M_op, eigvals(M_op))
            push!(evec_solutions_M_op, eigvecs(M_op))

        end
        t2 = time()

    elseif analytic_resonant==true && numeric_resonant==false

        t1 = time()
        for i in 1:length(eigenvalue_array)

            # Defining shorthand matrix: M_op
            global M_op = M_op_analytic_resonant(eigenvalue_array[i],τ,kyrho,η_i,η_e,omn,β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating λ-eigenvalues and Φ eigenvectors
            push!(eval_solutions_M_op, eigvals(M_op))
            push!(evec_solutions_M_op, eigvecs(M_op))

        end
        t2 = time()

    elseif numeric_resonant==true && analytic_resonant==false

        t1 = time()
        for i in 1:length(eigenvalue_array)

            # Defining shorthand matrix: M_op
            global M_op = M_op_numeric_resonant(eigenvalue_array[i],τ,kyrho,η_i,η_e,omn,β)

            # Applying boundary conditions to MΦ=0
            # Setting Φ=0 at edges of the domain
            M_op[1,:] = zeros(nz0+1) # Zeroing out first row of M_op
            M_op[end,:] = zeros(nz0+1) # Zeroing out last row of M_op
            M_op[1,1] = 1.0 # 1*Φ_LB=0 (LB: left boundary)
            M_op[end,end] = 1.0 # 1*Φ_RB=0 (RB: right boundary)

            # Calculating λ-eigenvalues and Φ eigenvectors
            push!(eval_solutions_M_op, eigvals(M_op))
            push!(evec_solutions_M_op, eigvecs(M_op))

        end
        t2 = time()

    end

    if quiet_time==false
        print("\nTime elapsed for calculating eigenvalues and eigenvectors:\n", round(t2-t1, digits=2), " seconds\n")
    end

end

# Calculating and plotting KBM eigenfunction candidates after finding eigenvalues from condition-number functions above 

function M_op_efn_array_plot(; 
        compute_all_evecs=true, # computing eigenvectors for all computed eigenvalues
        plot_soln_index=1, # determines which solution is plotted, but all solutions are evaluated if above argument is true 
        single_soln_index=1, # computing and plotting eigenvector for single specific eigenvalue solution; if true, overrides above plotting argument
        row_min_max=[1 2], 
        col_min_max=[1 3], 
        linewidth=2, 
        linewidth_axes=2,
        norm=false,
        calculate_eigenvectors=true, # set to false only after eigenvectors have been calculated; useful for replotting results  
        plot=true,
        start_λ_index=nz0-7,
        auto_index=true)

    if compute_all_evecs==false
        plot_soln_index=single_soln_index
    end

    global soln_index_M_op = plot_soln_index

    if calculate_eigenvectors==true
        eval_evec_soln_condition_M_op(compute_all_evecs=compute_all_evecs, single_soln_index=single_soln_index)
    end

    if compute_all_evecs==true
        global min_λ_index = indexin(minimum(abs.(eval_solutions_M_op[plot_soln_index])), abs.(eval_solutions_M_op[plot_soln_index]))[1]
    else
        global min_λ_index = indexin(minimum(abs.(eval_solutions_M_op[1])), abs.(eval_solutions_M_op[1]))[1]
    end

    if auto_index==true
        start_λ_index = min_λ_index - (row_min_max[2]*col_min_max[2]-1)
    end

    print("\nλ-min index = ", min_λ_index, "\n")

    if gene_units_condition_M_op==true
        print("\nEigenvalue = ", round(eval_condition_M_op_gene_units[plot_soln_index], digits=3), "\n")
    else
        print("\nEigenvalue = ", round(eval_condition_M_op[plot_soln_index], digits=3), "\n")
    end
    beta_plot = beta_condition_M_op

    if plot==true

        f = Figure(size=(1600,1000))

        fontsize_theme = Theme(fontsize = 25)
        set_theme!(fontsize_theme)

        for i in row_min_max[1]:row_min_max[2]
            
            for j in col_min_max[1]:col_min_max[2]

                global ax = Axis(f[i, j],
                                 title = latexstring(L"Φ \; - \; β = ", round(100*beta_plot, digits=2),L"% \; - \; \textrm{#}",start_λ_index + (i-1)*col_min_max[2] + (j-1)), 
                                  xlabel = L"θ / π",
                                  ylabel = L"a.u.",  
                                  titlesize = 20, 
                                  xlabelsize = 30, 
                                  ylabelsize = 30,
                                  xticklabelsize = 15,
                                  yticklabelsize = 15,
                                  xgridvisible = :false, 
                                  ygridvisible = :false) 

                if compute_all_evecs==true
                    if norm == true
                        p = lines!(thp, abs.(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(abs.(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :black)
                        p = lines!(thp, real(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(real(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :red)
                        p = lines!(thp, imag(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(imag(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :blue)
                    elseif norm == false
                        p = lines!(thp, abs.(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :black)
                        p = lines!(thp, real(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :red)
                        p = lines!(thp, imag(evec_solutions_M_op[plot_soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :blue)
                    end
                else 
                    if norm == true
                        p = lines!(thp, abs.(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(abs.(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :black)
                        p = lines!(thp, real(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(real(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :red)
                        p = lines!(thp, imag(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(imag(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :blue)
                    elseif norm == false
                        p = lines!(thp, abs.(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :black)
                        p = lines!(thp, real(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :red)
                        p = lines!(thp, imag(evec_solutions_M_op[1][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :blue)
                    end
                end

                p = hlines!(0, linewidth = linewidth_axes, color = :black)
                p = vlines!(0, linewidth = linewidth_axes, color = :black)

                hlines!(0, color = :black)

            end

        end

        f

    end

end

# Calculating and plotting KBM eigenfunction candidates based on eigenvalue and parameter inputs provided by user 

function M_op_efn_array_plot_any_input(eigenvalue; 
        β=0.01, 
        kyrho_gene=0.2,
        omti=3.5, 
        omte=0,
        omn=1,
        τ=1,
        soln_index=1, 
        analytic_resonant=0,
        numeric_resonant=0,
        row_min_max=[1 2], 
        col_min_max=[1 3], 
        linewidth=2, 
        linewidth_axes=2,
        norm=false,
        calculate_eigenvectors=true,
        plot=true,
        start_λ_index=nz0-7,
        auto_index=true)

    if calculate_eigenvectors==true
        eval_evec_soln_M_op_any_input(eigenvalue,
                                      β=β, 
                                      kyrho_gene=kyrho_gene,
                                      omti=omti, 
                                      omte=omte,
                                      omn=omn,
                                      τ=τ,
                                      analytic_resonant=analytic_resonant,
                                      numeric_resonant=numeric_resonant)
    end

    global eigenvalue_M_op_any_input = eigenvalue
    global beta_M_op_any_input = β
    global soln_index_M_op_any_input = soln_index

    global min_λ_index = indexin(minimum(abs.(eval_solutions_M_op[soln_index])), abs.(eval_solutions_M_op[soln_index]))[1]

    if auto_index==true
        start_λ_index = min_λ_index - (row_min_max[2]*col_min_max[2]-1)
    end
            
    print("\nEigenvalue = ", round(eigenvalue[soln_index], digits=3), "\n")

    print("\nλ-min index = ", min_λ_index, "\n")

    if plot==true

        f = Figure(size=(1600,1000))

        fontsize_theme = Theme(fontsize = 25)
        set_theme!(fontsize_theme)

        for i in row_min_max[1]:row_min_max[2]
            
            for j in col_min_max[1]:col_min_max[2]

                global ax = Axis(f[i, j],
                                 title = latexstring(L"Φ \; - \; β = ", round(100*β, digits=2),L"% \; - \; \textrm{#}",start_λ_index + (i-1)*col_min_max[2] + (j-1)), 
                                  xlabel = L"θ / π",
                                  ylabel = L"a.u.",  
                                  titlesize = 20, 
                                  xlabelsize = 30, 
                                  ylabelsize = 30,
                                  xticklabelsize = 15,
                                  yticklabelsize = 15,
                                  xgridvisible = :false, 
                                  ygridvisible = :false) 

                if norm == true
                    p = lines!(thp, abs.(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(abs.(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :black)
                    p = lines!(thp, real(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(real(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :red)
                    p = lines!(thp, imag(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]) / maximum(imag(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)])), linewidth = linewidth, color = :blue)
                elseif norm == false
                    p = lines!(thp, abs.(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :black)
                    p = lines!(thp, real(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :red)
                    p = lines!(thp, imag(evec_solutions_M_op[soln_index][:, start_λ_index + (i-1)*col_min_max[2] + (j-1)]), linewidth = linewidth, color = :blue)
                end

                p = hlines!(0, linewidth = linewidth_axes, color = :black)
                p = vlines!(0, linewidth = linewidth_axes, color = :black)

                hlines!(0, color = :black)

            end

        end

        f

    end

end

function reset_single_M_op_efn_save()
    global single_M_op_efn_save_beta = Float64[]
    global single_M_op_efn_save_evals = ComplexF64[]
    global single_M_op_efn_save_full_evecs = Complex.(zeros(nz0+1))
end

# Plotting single eigenvector solution - can be run after using M_op_efn_array_plot() 

function single_M_op_efn_plot(
        λ_index; 
        soln_index=soln_index_M_op,
        plot=true,
        linewidth=2, 
        linewidth_axes=2,
        norm=false)

    if gene_units_condition_M_op==true
        print("\nEigenvalue = ", round(eval_condition_M_op_gene_units[soln_index], digits=3), "\n")
        append!(single_M_op_efn_save_evals, eval_condition_M_op_gene_units[soln_index])
    else
        print("\nEigenvalue = ", round(eval_condition_M_op[soln_index], digits=3), "\n")
        append!(single_M_op_efn_save_evals, eval_condition_M_op[soln_index])
    end
    beta_single_M_op=beta_condition_M_op

    if single_M_op_efn_save_full_evecs==Complex.(zeros(nz0+1))
        single_M_op_efn_save_full_evecs .= evec_solutions_M_op[soln_index][:, λ_index]
    else
        global single_M_op_efn_save_full_evecs = hcat(single_M_op_efn_save_full_evecs, evec_solutions_M_op[soln_index][:, λ_index])
    end

    append!(single_M_op_efn_save_beta, beta_single_M_op) 

    if plot==true

        f = Figure(size=(1600,1000))

        fontsize_theme = Theme(fontsize = 25)
        set_theme!(fontsize_theme)

        ax1  = Axis(f[1,1],
                                  
                    title = latexstring(L"Φ \; - \; β = ", round(100*beta_single_M_op, digits=2),L"% \; - \; \textrm{#}", λ_index), 
                    xlabel = L"θ / π",
                    ylabel = L"a.u.",  
                    titlesize = 30, 
                    xlabelsize = 40, 
                    ylabelsize = 40,
                    xticklabelsize = 30,
                    yticklabelsize = 30,
                    xgridvisible = :false, 
                    ygridvisible = :false) 
        
        if norm==true
            p = lines!(thp, abs.(evec_solutions_M_op[soln_index][:, λ_index]) / maximum(abs.(evec_solutions_M_op[soln_index][:, λ_index])), linewidth = linewidth, color = :black)
            p = lines!(thp, real(evec_solutions_M_op[soln_index][:, λ_index]) / maximum(real(evec_solutions_M_op[soln_index][:, λ_index])), linewidth = linewidth, color = :red)
            p = lines!(thp, imag(evec_solutions_M_op[soln_index][:, λ_index]) / maximum(imag(evec_solutions_M_op[soln_index][:, λ_index])), linewidth = linewidth, color = :blue)
        elseif norm==false
            p = lines!(thp, abs.(evec_solutions_M_op[soln_index][:, λ_index]), linewidth = linewidth, color = :black)
            p = lines!(thp, real(evec_solutions_M_op[soln_index][:, λ_index]), linewidth = linewidth, color = :red)
            p = lines!(thp, imag(evec_solutions_M_op[soln_index][:, λ_index]), linewidth = linewidth, color = :blue)
        end

        p = hlines!(0, linewidth = linewidth_axes, color = :black)
        p = vlines!(0, linewidth = linewidth_axes, color = :black)

        hlines!(0, color = :black)

        f

    end 

end

# Plotting single eigenvector solution - can be run after using M_op_efn_array_plot_any_input() 

function single_M_op_efn_plot_any_input(
        λ_index; 
        soln_index=soln_index_M_op_any_input,
        plot=true,
        linewidth=2, 
        linewidth_axes=2,
        norm=false)

    eigenvalue_single_M_op = eigenvalue_M_op_any_input
    beta_single_M_op = beta_M_op_any_input

    print("\nEigenvalue = ", round(eigenvalue_single_M_op[soln_index], digits=3), "\n")

    # if statement checks to see if initialised array is empty or partially full
    if single_M_op_efn_save_full_evecs==Complex.(zeros(nz0+1))
        single_M_op_efn_save_full_evecs .= evec_solutions_M_op[soln_index][:, λ_index]
    else
        global single_M_op_efn_save_full_evecs = hcat(single_M_op_efn_save_full_evecs, evec_solutions_M_op[soln_index][:, λ_index])
    end

    append!(single_M_op_efn_save_beta, beta_single_M_op) 

    if plot==true

        f = Figure(size=(1600,1000))

        fontsize_theme = Theme(fontsize = 25)
        set_theme!(fontsize_theme)

        ax1  = Axis(f[1,1],
                                  
                    title = latexstring(L"Φ \; - \; β = ", round(100*beta_single_M_op, digits=2),L"% \; - \; \textrm{#}", λ_index), 
                    xlabel = L"θ / π",
                    ylabel = L"a.u.",  
                    titlesize = 30, 
                    xlabelsize = 40, 
                    ylabelsize = 40,
                    xticklabelsize = 30,
                    yticklabelsize = 30,
                    xgridvisible = :false, 
                    ygridvisible = :false) 
        
        if norm==true
            p = lines!(thp, abs.(evec_solutions_M_op[soln_index][:, λ_index]) / maximum(abs.(evec_solutions_M_op[soln_index][:, λ_index])), linewidth = linewidth, color = :black)
            p = lines!(thp, real(evec_solutions_M_op[soln_index][:, λ_index]) / maximum(real(evec_solutions_M_op[soln_index][:, λ_index])), linewidth = linewidth, color = :red)
            p = lines!(thp, imag(evec_solutions_M_op[soln_index][:, λ_index]) / maximum(imag(evec_solutions_M_op[soln_index][:, λ_index])), linewidth = linewidth, color = :blue)
        elseif norm==false
            p = lines!(thp, abs.(evec_solutions_M_op[soln_index][:, λ_index]), linewidth = linewidth, color = :black)
            p = lines!(thp, real(evec_solutions_M_op[soln_index][:, λ_index]), linewidth = linewidth, color = :red)
            p = lines!(thp, imag(evec_solutions_M_op[soln_index][:, λ_index]), linewidth = linewidth, color = :blue)
        end

        p = hlines!(0, linewidth = linewidth_axes, color = :black)
        p = vlines!(0, linewidth = linewidth_axes, color = :black)

        hlines!(0, color = :black)

        f

    end 

end

##################################### End ######################################
