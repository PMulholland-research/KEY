############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Finding single eigenvalue, always starting with wide search-mesh and zooms until converged (repeats for each β)
# Use this for GIST inputs

function automatic_eigenvalue_finder_general_geometry(;
        beta=0.01,
        omti=3.5,
        omte=0,
        omn=1,
        kyrho_gene=0.2,
        n_ev=1, # number of eigenvalues to save
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_start=[0,1], # providing starting value for eigenvalue search
        omega_lims_start=[0,1], # providing starting value for eigenvalue search
        gamma_step_start=0.1,
        omega_step_start=0.1,
        grid_factor=0.1, # smaller number → finer grid in subsequent search box
        plot=false,
        import_geo=true, # set to false if one wishes to use already imported geometry (for inconsistent betascans)
        geo_dir="./input/",
        geo_name="filename.dat",
        gyy_rescale=false,
        gyy_rescale_factor=1,
        Ky_rescale=false,
        Ky_rescale_factor=1)

    global gamma_step=gamma_step_start
    global omega_step=omega_step_start

    global gamma_lims=gamma_lims_start
    global omega_lims=omega_lims_start

    # Initialising matrices for saving output
    global eigenvalue_save_converge=ComplexF64[100000+100000im, 10000+10000im, 1000+1000im, 100+100im, 10+10im]
        
    if import_geo==true
        import_gist(string(geo_dir,geo_name),
                    gyy_rescale=gyy_rescale,
                    gyy_rescale_factor=gyy_rescale_factor,
                    Ky_rescale=Ky_rescale,
                    Ky_rescale_factor=Ky_rescale_factor)
    end

    while abs(eigenvalue_save_converge[end]-eigenvalue_save_converge[end-1]) > 1e-10

        condition_M_op(beta,
                       gamma_lims=gamma_lims,
                       omega_lims=omega_lims,
                       gamma_step=gamma_step, 
                       omega_step=omega_step, 
                       kyrho_gene=kyrho_gene,
                       omti=omti,
                       omte=omte,
                       omn=omn,
                       n_ev=n_ev,
                       plot=plot,
                       analytic_resonant=analytic_resonant,
                       numeric_resonant=numeric_resonant,
                       gene_units=1)

        append!(eigenvalue_save_converge, eval_condition_M_op_gene_units)

        # Resetting starting point to be in vicinity of last-found eigenvalue

        global gamma_restart=imag(eval_condition_M_op_gene_units[1])
        global omega_restart=real(eval_condition_M_op_gene_units[1])

        if abs(gamma_restart)==0.0 && abs(omega_restart)==0.0
            global gamma_lims=[0,0.1]
            global omega_lims=[0,0.1]
        elseif abs(omega_restart)==0.0 
            global gamma_lims=[gamma_restart-gamma_step,gamma_restart+gamma_step]
            global omega_lims=[0,0.1]
        elseif abs(gamma_restart)==0.0 
            global gamma_lims=[0,0.1]
            global omega_lims=[omega_restart-omega_step,omega_restart+omega_step]
        else
           global gamma_lims=[gamma_restart-gamma_step,gamma_restart+gamma_step]
           global omega_lims=[omega_restart-omega_step,omega_restart+omega_step]
        end

        global gamma_step = grid_factor*(gamma_lims[2]-gamma_lims[1])
        global omega_step = grid_factor*(omega_lims[2]-omega_lims[1])

        # Failsafe to break out of for loop if solver bounces back and forth between two solutions
        if eigenvalue_save_converge[end]==eigenvalue_save_converge[end-2] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        elseif eigenvalue_save_converge[end]==eigenvalue_save_converge[end-3] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        elseif eigenvalue_save_converge[end]==eigenvalue_save_converge[end-4] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        elseif eigenvalue_save_converge[end]==eigenvalue_save_converge[end-5] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        end

    end

end

# Looping over eigenvalue finder (zooming until converged on solution) 
# Best option for dominant eigenvalues/eigenvectors

function automatic_betascan_general_geometry_zoomer(;
        omti=3.5,
        omte=0,
        omn=1,
        kyrho_gene=0.2,
        npol=20,
        gridpoints_per_pol=128,
        n_ev=1, # number of eigenvalues to find
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_start=[0,2], # providing starting value for eigenvalue search
        omega_lims_start=[0,1], # providing starting value for eigenvalue search
        gamma_step_start=0.1,
        omega_step_start=0.1,
        grid_factor=0.1, # smaller number → finer grid in subsequent search box
        eigenvector=false,
        import_geo=true, # set to false if one wishes to use already imported geometry (for inconsistent betascans)
        geo_dir="./input/",
        geo_name="filename.dat",
        padding=3,
        gyy_rescale=false,
        gyy_rescale_factor=1,
        Ky_rescale=false,
        Ky_rescale_factor=1,
        beta_from_file=true,
        beta_data_file="beta.dat",
        beta_label_file="beta_label.dat",
        beta_index=1:1:14,
        beta_input=0.01:0.01:0.05,
        save=false,
        savefolder="./output/")

    if analytic_resonant==false && numeric_resonant==false
        model="non_resonant"
    elseif analytic_resonant==true
        model="analytic_resonant"
    elseif numeric_resonant==true
        model="numeric_resonant"
    end

    global nz0=gridpoints_per_pol*npol

    if beta_from_file==true
        global all_beta=vec((readdlm(string(geo_dir,beta_data_file))))
        global all_beta_label=lpad.(vec(Int64.(readdlm(string(geo_dir,beta_label_file)))), padding, '0')
        global beta = all_beta[beta_index]
        global beta_label = all_beta_label[beta_index]
    else
        global beta = beta_input
        global beta_label = string.(beta)
    end

    if save==true
        # Creating directory for storing output data
        betascan_directory = string(savefolder,model,"_kyrho_gene_",kyrho_gene,"_betascan_",first(round.(beta,digits=4)),"_",last(round.(beta,digits=4)),"_n_ev_",n_ev,"/")
        mkdir(betascan_directory)
    end
    
    # Creating parameters.txt file for reproducibility
    global parameters=String[]
    push!(parameters, "# KEY parameters file")
    push!(parameters, "")
    push!(parameters, "# Geometry")
    push!(parameters, string("import_geo=",import_geo))
    push!(parameters, string("geo_directory=",geo_dir))
    push!(parameters, string("geo_name=",geo_name))
    push!(parameters, string("gyy_rescale=",gyy_rescale))
    push!(parameters, string("gyy_rescale_factor=",gyy_rescale_factor))
    push!(parameters, string("Ky_rescale=",Ky_rescale))
    push!(parameters, string("Ky_rescale_factor=",Ky_rescale_factor))
    push!(parameters, "")
    push!(parameters, "# Plasma")
    push!(parameters, string("beta=",round.(beta,digits=5)))
    push!(parameters, string("kyrho_gene=",kyrho_gene))
    push!(parameters, string("omti=",omti))
    push!(parameters, string("omte=",omte))
    push!(parameters, string("omn=",omn))
    push!(parameters, "")
    push!(parameters, "# Simulation")
    push!(parameters, string("n_ev=",n_ev))
    push!(parameters, string("analytic_resonant=",analytic_resonant))
    push!(parameters, string("numeric_resonant=",numeric_resonant))
    push!(parameters, string("eigenvector=",eigenvector))
    push!(parameters, string("gamma_lims_start=",gamma_lims_start))
    push!(parameters, string("omega_lims_start=",omega_lims_start))
    push!(parameters, string("gamma_step_start=",gamma_step_start))
    push!(parameters, string("omega_step_start=",omega_step_start))
    push!(parameters, string("grid_factor=",grid_factor))
    push!(parameters, "")

    if save==true
        writedlm(string(betascan_directory,"parameters.txt"), parameters)
    end
 
    # Initialising matrices for saving output
    global eigenvalue_save_matrix=Complex.(zeros(length(beta), n_ev))
    global beta_save_matrix=zeros(length(beta), n_ev)

    if eigenvector==true
        # Initialising eigenvector-saving function 
        reset_single_M_op_efn_save() 
    end

    for i in ProgressBar(1:length(beta))

        automatic_eigenvalue_finder_general_geometry(;
                                    beta=beta[i],
                                    omti=omti,
                                    omte=omte,
                                    omn=omn,
                                    kyrho_gene=kyrho_gene,
                                    n_ev=1, # number of eigenvalues to save
                                    analytic_resonant=analytic_resonant,
                                    numeric_resonant=numeric_resonant,
                                    gamma_lims_start=gamma_lims_start, 
                                    omega_lims_start=omega_lims_start, 
                                    gamma_step_start=gamma_step_start,
                                    omega_step_start=omega_step_start,
                                    grid_factor=grid_factor,
                                    import_geo=import_geo,
                                    geo_dir=geo_dir,
                                    geo_name=string(geo_name,beta_label[i],".dat"),
                                    gyy_rescale=gyy_rescale,
                                    gyy_rescale_factor=gyy_rescale_factor,
                                    Ky_rescale=Ky_rescale,
                                    Ky_rescale_factor=Ky_rescale_factor)

        eigenvalue_save_matrix[i,1:n_ev].=last(eigenvalue_save_converge)
            
        beta_save_matrix[i,1:n_ev].=beta[i]
            
        if eigenvector==true
            M_op_efn_array_plot(plot=false)
            single_M_op_efn_plot(min_λ_index, plot=0) # saving eigenvector of dominant eigenmode 
        end

    end
    
    # Creating dataframe for dominant eigenvalues
    scanlog_matrix=hcat(beta_save_matrix[:,1],imag(eigenvalue_save_matrix[:,1]),real(eigenvalue_save_matrix[:,1]))
    scanlog_column_names=["beta","gamma","omega"]
    global scanlog_dataframe=DataFrame(scanlog_matrix, scanlog_column_names)
      
    if eigenvector==true

        global eigenvector_column_names= String[]
        map(x -> push!(eigenvector_column_names, string("eigenvector",x)), 1:1:length(beta))

        # Creating dataframe for dominant eigenvectors
        eigenvector_matrix=single_M_op_efn_save_full_evecs
        global eigenvector_dataframe=DataFrame(eigenvector_matrix, eigenvector_column_names)

    end

    if save==true
        CSV.write(string(betascan_directory,"eigenvalues.dat"),
                  scanlog_dataframe,
                  delim='\t')

        if eigenvector==true
            CSV.write(string(betascan_directory,"eigenvectors.dat"),
                      eigenvector_dataframe,
                      delim='\t')
        end 

    end

end

# Alternative to above: this betascan uses fixed search box to find eigenvalues 
# Possibly better option for n_ev > 1 (searching for subdominant modes) as it doesn't zoom in on dominant eigenvalue

function automatic_betascan_general_geometry_fixed_grid(;
        omti=3.5,
        omte=0,
        omn=1,
        kyrho_gene=0.2,
        npol=20,
        gridpoints_per_pol=128,
        n_ev=15, # number of eigenvalues to find
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_scan=[0,0.2], 
        omega_lims_scan=[0.4,1.2], 
        gamma_step_scan=0.01,
        omega_step_scan=0.01,
        eigenvector=false,
        import_geo=true, # set to false if one wishes to use already imported geometry (for inconsistent betascans)
        geo_dir="./input/",
        geo_name="filename.dat",
        padding=3,
        gyy_rescale=false,
        gyy_rescale_factor=1,
        Ky_rescale=false,
        Ky_rescale_factor=1,
        beta_from_file=true,
        beta_data_file="beta.dat",
        beta_label_file="beta_label.dat",
        beta_index=1:1:10,
        beta_input=0.01:0.01:0.05,
        save=false,
        savefolder="./output/")

    if analytic_resonant==false && numeric_resonant==false
        model="non_resonant"
    elseif analytic_resonant==true
        model="analytic_resonant"
    elseif numeric_resonant==true
        model="numeric_resonant"
    end

    global nz0=gridpoints_per_pol*npol

    if beta_from_file==true
        global all_beta=vec((readdlm(string(geo_dir,beta_data_file))))
        global all_beta_label=lpad.(vec(Int64.(readdlm(string(geo_dir,beta_label_file)))), padding, '0')
        global beta = vec(all_beta[beta_index])
        global beta_label = all_beta_label[beta_index]
    else
        global beta = beta_input
    end

    if save==true
        # Creating directory for storing output data
        betascan_directory = string(savefolder,model,"_kyrho_gene_",kyrho_gene,"_betascan_",first(round.(beta,digits=4)),"_",last(round.(beta,digits=4)),"_n_ev_",n_ev,"/")
        mkdir(betascan_directory)
    end
    
    # Creating parameters.txt file for reproducibility
    global parameters=String[]
    push!(parameters, "# KEY parameters file")
    push!(parameters, "")
    push!(parameters, "# Geometry")
    push!(parameters, string("import_geo=",import_geo))
    push!(parameters, string("geo_directory=",geo_dir))
    push!(parameters, string("geo_name=",geo_name))
    push!(parameters, string("gyy_rescale=",gyy_rescale))
    push!(parameters, string("gyy_rescale_factor=",gyy_rescale_factor))
    push!(parameters, string("Ky_rescale=",Ky_rescale))
    push!(parameters, string("Ky_rescale_factor=",Ky_rescale_factor))
    push!(parameters, "")
    push!(parameters, "# Plasma")
    push!(parameters, string("beta=",round.(beta,digits=4)))
    push!(parameters, string("kyrho_gene=",kyrho_gene))
    push!(parameters, string("omti=",omti))
    push!(parameters, string("omte=",omte))
    push!(parameters, string("omn=",omn))
    push!(parameters, "")
    push!(parameters, "# Simulation")
    push!(parameters, string("n_ev=",n_ev))
    push!(parameters, string("analytic_resonant=",analytic_resonant))
    push!(parameters, string("numeric_resonant=",numeric_resonant))
    push!(parameters, string("eigenvector=",eigenvector))
    push!(parameters, string("gamma_lims_scan=",gamma_lims_scan))
    push!(parameters, string("omega_lims_scan=",omega_lims_scan))
    push!(parameters, string("gamma_step_scan=",gamma_step_scan))
    push!(parameters, string("omega_step_scan=",omega_step_scan))
    push!(parameters, "")

    if save==true
        writedlm(string(betascan_directory,"parameters.txt"), parameters)
    end

    global gamma_step=gamma_step_scan
    global omega_step=omega_step_scan

    global gamma_lims=gamma_lims_scan
    global omega_lims=omega_lims_scan

    # Initialising matrices for saving output
    global eigenvalue_save_matrix=Complex.(zeros(length(beta), n_ev))
    global beta_save_matrix=zeros(length(beta), n_ev)

    if eigenvector==true
        # Initialising eigenvector-saving function 
        reset_single_M_op_efn_save() 
    end
        
    for i in ProgressBar(1:length(beta))

        if import_geo==true
            import_gist(string(geo_dir,geo_name,beta_label[i],".dat"),
                        gyy_rescale=gyy_rescale,
                        gyy_rescale_factor=gyy_rescale_factor,
                        Ky_rescale=Ky_rescale,
                        Ky_rescale_factor=Ky_rescale_factor)
        end

        condition_M_op(beta[i],
                       gamma_lims=gamma_lims,
                       omega_lims=omega_lims,
                       gamma_step=gamma_step, 
                       omega_step=omega_step, 
                       kyrho_gene=kyrho_gene,
                       omti=omti,
                       omte=omte,
                       omn=omn,
                       n_ev=n_ev,
                       plot=false,
                       analytic_resonant=analytic_resonant,
                       numeric_resonant=numeric_resonant,
                       gene_units=1)

        eigenvalue_save_matrix[i,1:n_ev]=eval_condition_M_op_gene_units
        beta_save_matrix[i,1:n_ev].=beta_condition_M_op 

    end
    
    # Creating dataframe for subdominant eigenspectrum: rows correspond to changing beta 
    eigenspectrum_matrix=eigenvalue_save_matrix
    eigenspectrum_column_names=[]
    map(x -> push!(eigenspectrum_column_names, string("EV",x)), 1:1:n_ev)

    global eigenspectrum_dataframe=DataFrame(eigenspectrum_matrix, eigenspectrum_column_names)

    if save==true
        CSV.write(string(betascan_directory,"eigenspectrum.dat"),
                  eigenspectrum_dataframe,
                  delim='\t')
    end

end

# Automatic ky-scan

function automatic_kyscan_general_geometry_zoomer(;
        kyrho_gene=0.05:0.05:0.8,
        omti=3.5,
        omte=0,
        omn=1,
        npol=20,
        gridpoints_per_pol=128,
        n_ev=1, # number of eigenvalues to find
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_start=[0,2], # providing starting value for eigenvalue search
        omega_lims_start=[0,1], # providing starting value for eigenvalue search
        gamma_step_start=0.1,
        omega_step_start=0.1,
        grid_factor=0.1, # smaller number → finer grid in subsequent search box
        eigenvector=false,
        import_geo=false, # set to false if one wishes to use already imported geometry (for inconsistent betascans)
        geo_dir="./input/",
        geo_name="filename.dat",
        padding=3,
        gyy_rescale=false,
        gyy_rescale_factor=1,
        Ky_rescale=false,
        Ky_rescale_factor=1,
        beta_from_file=false,
        beta_data_file="beta.dat",
        beta_label_file="beta_label.dat",
        beta_index=1,
        beta_input=0.05,
        save=false,
        savefolder="./output/")

    if analytic_resonant==false && numeric_resonant==false
        model="non_resonant"
    elseif analytic_resonant==true
        model="analytic_resonant"
    elseif numeric_resonant==true
        model="numeric_resonant"
    end

    global nz0=gridpoints_per_pol*npol

    if beta_from_file==true
        global all_beta=vec((readdlm(string(geo_dir,beta_data_file))))
        global all_beta_label=lpad.(vec(Int64.(readdlm(string(geo_dir,beta_label_file)))), padding, '0')
        global beta = all_beta[beta_index]
        global beta_label = all_beta_label[beta_index]
    else
        global beta = beta_input
        global beta_label = string.(beta)
    end

    if save==true
        # Creating directory for storing output data
        kyscan_directory = string(savefolder,model,"_beta_",beta,"_kyrho_gene_scan_",first(kyrho_gene),"_",last(kyrho_gene),"_n_ev_",n_ev,"/")
        mkdir(betascan_directory)
    end
    
    # Creating parameters.txt file for reproducibility
    global parameters=String[]
    push!(parameters, "# KEY parameters file")
    push!(parameters, "")
    push!(parameters, "# Geometry")
    push!(parameters, string("import_geo=",import_geo))
    push!(parameters, string("geo_directory=",geo_dir))
    push!(parameters, string("geo_name=",geo_name))
    push!(parameters, string("gyy_rescale=",gyy_rescale))
    push!(parameters, string("gyy_rescale_factor=",gyy_rescale_factor))
    push!(parameters, string("Ky_rescale=",Ky_rescale))
    push!(parameters, string("Ky_rescale_factor=",Ky_rescale_factor))
    push!(parameters, "")
    push!(parameters, "# Plasma")
    push!(parameters, string("beta=",beta))
    push!(parameters, string("kyrho_gene=",kyrho_gene))
    push!(parameters, string("omti=",omti))
    push!(parameters, string("omte=",omte))
    push!(parameters, string("omn=",omn))
    push!(parameters, "")
    push!(parameters, "# Simulation")
    push!(parameters, string("n_ev=",n_ev))
    push!(parameters, string("analytic_resonant=",analytic_resonant))
    push!(parameters, string("numeric_resonant=",numeric_resonant))
    push!(parameters, string("eigenvector=",eigenvector))
    push!(parameters, string("gamma_lims_start=",gamma_lims_start))
    push!(parameters, string("omega_lims_start=",omega_lims_start))
    push!(parameters, string("gamma_step_start=",gamma_step_start))
    push!(parameters, string("omega_step_start=",omega_step_start))
    push!(parameters, string("grid_factor=",grid_factor))
    push!(parameters, "")

    if save==true
        writedlm(string(betascan_directory,"parameters.txt"), parameters)
    end
 
    # Initialising matrices for saving output
    global eigenvalue_save_matrix=Complex.(zeros(length(kyrho_gene), n_ev))
    global ky_save_matrix=zeros(length(kyrho_gene), n_ev)

    if eigenvector==true
        # Initialising eigenvector-saving function 
        reset_single_M_op_efn_save() 
    end

    for i in ProgressBar(1:length(kyrho_gene))

        automatic_eigenvalue_finder_general_geometry(;
                                    beta=beta,
                                    omti=omti,
                                    omte=omte,
                                    omn=omn,
                                    kyrho_gene=kyrho_gene[i],
                                    n_ev=1, # number of eigenvalues to save
                                    analytic_resonant=analytic_resonant,
                                    numeric_resonant=numeric_resonant,
                                    gamma_lims_start=gamma_lims_start, 
                                    omega_lims_start=omega_lims_start, 
                                    gamma_step_start=gamma_step_start,
                                    omega_step_start=omega_step_start,
                                    grid_factor=grid_factor,
                                    import_geo=import_geo,
                                    geo_dir=geo_dir,
                                    geo_name=string(geo_name,beta_label,".dat"),
                                    gyy_rescale=gyy_rescale,
                                    gyy_rescale_factor=gyy_rescale_factor,
                                    Ky_rescale=Ky_rescale,
                                    Ky_rescale_factor=Ky_rescale_factor)


        eigenvalue_save_matrix[i,1:n_ev].=last(eigenvalue_save_converge)
            
        ky_save_matrix[i,1:n_ev].=kyrho_gene[i]
            
        if eigenvector==true
            M_op_efn_array_plot(plot=false)
            single_M_op_efn_plot(min_λ_index, plot=0) # saving eigenvector of dominant eigenmode 
        end

    end
    
    # Creating dataframe for dominant eigenvalues
    scanlog_matrix=hcat(ky_save_matrix[:,1],imag(eigenvalue_save_matrix[:,1]),real(eigenvalue_save_matrix[:,1]))
    scanlog_column_names=["kyrho_gene","gamma","omega"]
    global scanlog_dataframe=DataFrame(scanlog_matrix, scanlog_column_names)
      
    if eigenvector==true

        global eigenvector_column_names= String[]
        #map(x -> push!(eigenvector_labels, string("Φ (β=",100*beta[x],"%)")), 1:1:length(beta))
        map(x -> push!(eigenvector_column_names, string("eigenvector",x)), 1:1:length(beta))

        # Creating dataframe for dominant eigenvectors
        eigenvector_matrix=single_M_op_efn_save_full_evecs
        global eigenvector_dataframe=DataFrame(eigenvector_matrix, eigenvector_column_names)

    end

    if save==true
        CSV.write(string(kyscan_directory,"eigenvalues.dat"),
                  scanlog_dataframe,
                  delim='\t')

        if eigenvector==true
            CSV.write(string(kyscan_directory,"eigenvectors.dat"),
                      eigenvector_dataframe,
                      delim='\t')
        end 

    end

end

# The following functions are the same as above, but now for s-alpha
# Automating the calculation of eigenvalues/eigenvectors for β/k_y scan in s-alpha

# Finding single eigenvalue, always starting with wide search-mesh and zooms until convereged (repeats for each β)
function automatic_eigenvalue_finder_salpha(;
        beta=0.01,
        shear=0.796,
        npol=20,
        gridpoints_per_pol=64,
        omti=10,
        omte=10,
        omn=10,
        import_geo=true, # set to false if one wishes to use already imported geometry (for inconsistent betascans)
        kyrho_gene=0.1,
        n_ev=1, # number of eigenvalues to save
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_start=[0,6], # providing starting value for eigenvalue search
        omega_lims_start=[0,6], # providing starting value for eigenvalue search
        gamma_step_start=0.1,
        omega_step_start=0.1,
        grid_factor=0.1, # smaller number → finer grid in subsequent search box
        progress_bar=true,
        quiet_output=false,
        plot=false)

    global gamma_step=gamma_step_start
    global omega_step=omega_step_start

    global gamma_lims=gamma_lims_start
    global omega_lims=omega_lims_start

    # Initialising matrices for saving output
    global eigenvalue_save_converge=ComplexF64[100000+100000im, 10000+10000im, 1000+1000im, 100+100im, 10+10im]
        
    if import_geo==true

    global nz0=gridpoints_per_pol*npol

        import_salpha(shear=shear,
                      npol=npol,
                      gridpoints=nz0,
                      beta=beta,
                      omti=omti,
                      omte=omte,
                      omn=omn)
    end

    while abs(eigenvalue_save_converge[end]-eigenvalue_save_converge[end-1]) > 1e-10

        condition_M_op(beta,
                       gamma_lims=gamma_lims,
                       omega_lims=omega_lims,
                       gamma_step=gamma_step, 
                       omega_step=omega_step, 
                       kyrho_gene=kyrho_gene,
                       omti=omti,
                       omte=omte,
                       omn=omn,
                       n_ev=n_ev,
                       plot=plot,
                       analytic_resonant=analytic_resonant,
                       numeric_resonant=numeric_resonant,
                       gene_units=1,
                       progress_bar=progress_bar,
                       quiet_output=quiet_output)

        append!(eigenvalue_save_converge, eval_condition_M_op_gene_units)

        # Resetting starting point to be in vicinity of last-found eigenvalue
        
        global gamma_restart=imag(eval_condition_M_op_gene_units[1])
        global omega_restart=real(eval_condition_M_op_gene_units[1])

        if abs(gamma_restart)==0.0 && abs(omega_restart)==0.0
            global gamma_lims=[0,0.1]
            global omega_lims=[0,0.1]
        elseif abs(omega_restart)==0.0 
            global gamma_lims=[gamma_restart-gamma_step,gamma_restart+gamma_step]
            global omega_lims=[0,0.1]
        elseif abs(gamma_restart)==0.0 
            global gamma_lims=[0,0.1]
            global omega_lims=[omega_restart-omega_step,omega_restart+omega_step]
        else
           global gamma_lims=[gamma_restart-gamma_step,gamma_restart+gamma_step]
           global omega_lims=[omega_restart-omega_step,omega_restart+omega_step]
        end

        global gamma_step = grid_factor*(gamma_lims[2]-gamma_lims[1])
        global omega_step = grid_factor*(omega_lims[2]-omega_lims[1])

        # Failsafe to break out of for loop if solver bounces back and forth between two solutions
        if eigenvalue_save_converge[end]==eigenvalue_save_converge[end-2] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        elseif eigenvalue_save_converge[end]==eigenvalue_save_converge[end-3] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        elseif eigenvalue_save_converge[end]==eigenvalue_save_converge[end-4] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        elseif eigenvalue_save_converge[end]==eigenvalue_save_converge[end-5] 
            print("\nYou're stuck in a loop! Try adjusting your search parameters.\n")
            break
        end

    end

end

# Looping over eigenvalue finder (zooming until converged on solution) 
# Best option for dominant eigenvalues/eigenvectors

function automatic_betascan_salpha_zoomer(;
        betascan=0.001:0.001:0.05,
        shear=0.796,
        npol=20,
        gridpoints_per_pol=64,
        omti=10,
        omte=10,
        omn=10,
        import_geo=true, # set to false if one wishes to use already imported geometry (for inconsistent betascans)
        kyrho_gene=0.1,
        n_ev=1, # number of eigenvalues to find
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_start=[0,7], # providing starting value for eigenvalue search
        omega_lims_start=[0,4], # providing starting value for eigenvalue search
        gamma_step_start=0.1,
        omega_step_start=0.1,
        grid_factor=0.1, # smaller number → finer grid in subsequent search box
        eigenvector=false,
        progress_bar=true,
        quiet_output=false,
        quiet_scan=false,
        save=false,
        savefolder="./output/")

    if analytic_resonant==false && numeric_resonant==false
        model="non_resonant"
    elseif analytic_resonant==true
        model="analytic_resonant"
    elseif numeric_resonant==true
        model="numeric_resonant"
    end

    global nz0=gridpoints_per_pol*npol

    global beta=betascan

    if save==true
        # Creating directory for storing output data
        betascan_directory = string(savefolder,"salpha_shear_",shear,"_",model,"_kyrho_gene_",kyrho_gene,"_betascan_",first(round.(beta,digits=4)),"_",last(round.(beta,digits=4)),"_n_ev_",n_ev,"/")
        mkdir(betascan_directory)
    end
    
    # Creating parameters.txt file for reproducibility
    global parameters=String[]
    push!(parameters, "# KEY parameters file")
    push!(parameters, "")
    push!(parameters, "# Geometry")
    push!(parameters, "geo_name=salpha")
    push!(parameters, string("import_geo=",import_geo))
    push!(parameters, string("shear=",shear))
    push!(parameters, string("npol=",npol))
    push!(parameters, string("gridpoints=",nz0))
    push!(parameters, "")
    push!(parameters, "# Plasma")
    push!(parameters, string("beta=",round.(beta,digits=5)))
    push!(parameters, string("kyrho_gene=",kyrho_gene))
    push!(parameters, string("omti=",omti))
    push!(parameters, string("omte=",omte))
    push!(parameters, string("omn=",omn))
    push!(parameters, "")
    push!(parameters, "# Simulation")
    push!(parameters, string("n_ev=",n_ev))
    push!(parameters, string("analytic_resonant=",analytic_resonant))
    push!(parameters, string("numeric_resonant=",numeric_resonant))
    push!(parameters, string("eigenvector=",eigenvector))
    push!(parameters, string("gamma_lims_start=",gamma_lims_start))
    push!(parameters, string("omega_lims_start=",omega_lims_start))
    push!(parameters, string("gamma_step_start=",gamma_step_start))
    push!(parameters, string("omega_step_start=",omega_step_start))
    push!(parameters, string("grid_factor=",grid_factor))
    push!(parameters, "")

    if save==true
        writedlm(string(betascan_directory,"parameters.txt"), parameters)
    end
 
    # Initialising matrices for saving output
    global eigenvalue_save_matrix=Complex.(zeros(length(beta), n_ev))
    global beta_save_matrix=zeros(length(beta), n_ev)

    if eigenvector==true
        # Initialising eigenvector-saving function 
        reset_single_M_op_efn_save() 
    end
    
    # Quiet and loud options: with and without progress bar
    progress_bar_array = ProgressBar(1:length(beta)) 
    for_loop_array = 1:length(beta) 

    if quiet_scan==false
        scan_array = progress_bar_array
    else
        scan_array = for_loop_array
    end

    for i in scan_array

        automatic_eigenvalue_finder_salpha(;
                                           beta=beta[i],
                                           shear=shear,
                                           npol=npol,
                                           gridpoints_per_pol=gridpoints_per_pol,
                                           omti=omti,
                                           omte=omte,
                                           omn=omn,
                                           import_geo=import_geo,
                                           kyrho_gene=kyrho_gene,
                                           n_ev=1, # number of eigenvalues to save
                                           analytic_resonant=analytic_resonant,
                                           numeric_resonant=numeric_resonant,
                                           gamma_lims_start=gamma_lims_start, 
                                           omega_lims_start=omega_lims_start, 
                                           gamma_step_start=gamma_step_start,
                                           omega_step_start=omega_step_start,
                                           grid_factor=grid_factor,
                                           progress_bar=progress_bar,
                                           quiet_output=quiet_output)

        eigenvalue_save_matrix[i,1:n_ev].=last(eigenvalue_save_converge)
            
        beta_save_matrix[i,1:n_ev].=beta[i]
            
        if eigenvector==true
            M_op_efn_array_plot(plot=false)
            single_M_op_efn_plot(min_λ_index, plot=0) # saving eigenvector of dominant eigenmode 
        end

    end
    
    # Creating dataframe for dominant eigenvalues
    scanlog_matrix=hcat(beta_save_matrix[:,1],imag(eigenvalue_save_matrix[:,1]),real(eigenvalue_save_matrix[:,1]))
    scanlog_column_names=["beta","gamma","omega"]
    global scanlog_dataframe=DataFrame(scanlog_matrix, scanlog_column_names)
      
    if eigenvector==true

        global eigenvector_column_names= String[]
        map(x -> push!(eigenvector_column_names, string("eigenvector",x)), 1:1:length(beta))

        # Creating dataframe for dominant eigenvectors
        eigenvector_matrix=single_M_op_efn_save_full_evecs
        global eigenvector_dataframe=DataFrame(eigenvector_matrix, eigenvector_column_names)

    end

    if save==true
        CSV.write(string(betascan_directory,"eigenvalues.dat"),
                  scanlog_dataframe,
                  delim='\t')

        if eigenvector==true
            CSV.write(string(betascan_directory,"eigenvectors.dat"),
                      eigenvector_dataframe,
                      delim='\t')
        end 

    end

end

# Alternative to above: this betascan uses fixed search box to find eigenvalues 
# Possibly better option for n_ev > 1 (searching for subdominant modes) as it doesn't zoom in on dominant eigenvalue

function automatic_betascan_salpha_fixed_grid(;
        beta=0.01:0.01:0.05,
        shear=0.8,
        npol=20,
        gridpoints_per_pol=64,
        omti=10,
        omte=10,
        omn=10,
        import_geo=true,
        kyrho_gene=0.1,
        n_ev=10, # number of eigenvalues to find
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_scan=[0,1], 
        omega_lims_scan=[0,1], 
        gamma_step_scan=0.1,
        omega_step_scan=0.1,
        eigenvector=false,
        save=true,
        savefolder="./output/")

    if analytic_resonant==false && numeric_resonant==false
        model="non_resonant"
    elseif analytic_resonant==true
        model="analytic_resonant"
    elseif numeric_resonant==true
        model="numeric_resonant"
    end

    global nz0=gridpoints_per_pol*npol

    global beta=beta

    if save==true
        # Creating directory for storing output data
        betascan_directory = string(savefolder,"salpha_shear_",shear,"_",model,"_kyrho_gene_",kyrho_gene,"_betascan_",first(round.(beta,digits=4)),"_",last(round.(beta,digits=4)),"_n_ev_",n_ev,"/")
        mkdir(betascan_directory)
    end
    
    # Creating parameters.txt file for reproducibility
    global parameters=String[]
    push!(parameters, "# KEY parameters file")
    push!(parameters, "")
    push!(parameters, "# Geometry")
    push!(parameters, "geo_name=salpha")
    push!(parameters, string("import_geo=",import_geo))
    push!(parameters, string("shear=",shear))
    push!(parameters, string("npol=",npol))
    push!(parameters, string("gridpoints=",nz0))
    push!(parameters, "")
    push!(parameters, "# Plasma")
    push!(parameters, string("beta=",round.(beta,digits=5)))
    push!(parameters, string("kyrho_gene=",kyrho_gene))
    push!(parameters, string("omti=",omti))
    push!(parameters, string("omte=",omte))
    push!(parameters, string("omn=",omn))
    push!(parameters, "")
    push!(parameters, "# Simulation")
    push!(parameters, string("n_ev=",n_ev))
    push!(parameters, string("analytic_resonant=",analytic_resonant))
    push!(parameters, string("numeric_resonant=",numeric_resonant))
    push!(parameters, string("eigenvector=",eigenvector))
    push!(parameters, string("gamma_lims_scan=",gamma_lims_scan))
    push!(parameters, string("omega_lims_scan=",omega_lims_scan))
    push!(parameters, string("gamma_step_scan=",gamma_step_scan))
    push!(parameters, string("omega_step_scan=",omega_step_scan))
    push!(parameters, "")

    if save==true
        writedlm(string(betascan_directory,"parameters.txt"), parameters)
    end

    global gamma_step=gamma_step_scan
    global omega_step=omega_step_scan

    global gamma_lims=gamma_lims_scan
    global omega_lims=omega_lims_scan

    # Initialising matrices for saving output
    global eigenvalue_save_matrix=Complex.(zeros(length(beta), n_ev))
    global beta_save_matrix=zeros(length(beta), n_ev)

    if eigenvector==true
        # Initialising eigenvector-saving function 
        reset_single_M_op_efn_save() 
    end
        
    for i in ProgressBar(1:length(beta))
    
        if import_geo==true

            import_salpha(shear=shear,
                          npol=npol,
                          gridpoints=nz0,
                          beta=beta[i],
                          omti=omti,
                          omte=omte,
                          omn=omn)
        end

        condition_M_op(beta[i],
                       gamma_lims=gamma_lims,
                       omega_lims=omega_lims,
                       gamma_step=gamma_step, 
                       omega_step=omega_step, 
                       kyrho_gene=kyrho_gene,
                       omti=omti,
                       omte=omte,
                       omn=omn,
                       n_ev=n_ev,
                       plot=false,
                       analytic_resonant=analytic_resonant,
                       numeric_resonant=numeric_resonant,
                       gene_units=1)

        eigenvalue_save_matrix[i,1:n_ev]=eval_condition_M_op_gene_units
        beta_save_matrix[i,1:n_ev].=beta_condition_M_op 

    end

    # Creating dataframe for subdominant eigenspectrum: rows correspond to changing beta 
    eigenspectrum_matrix=eigenvalue_save_matrix
    eigenspectrum_column_names=[]
    map(x -> push!(eigenspectrum_column_names, string("EV",x)), 1:1:n_ev)

    global eigenspectrum_dataframe=DataFrame(eigenspectrum_matrix, eigenspectrum_column_names)

    if save==true
        CSV.write(string(betascan_directory,"eigenspectrum.dat"),
                  eigenspectrum_dataframe,
                  delim='\t')
    end

end

# Automatic ky-scan

function automatic_kyscan_salpha_zoomer(;
        kyrho_gene=0.05:0.05:0.8,
        beta=0.01,
        shear=0.8,
        npol=20,
        gridpoints_per_pol=64,
        omti=10,
        omte=10,
        omn=10,
        import_geo=true, # set to false if one wishes to use already imported geometry (for inconsistent betascans)
        n_ev=1, # number of eigenvalues to find
        analytic_resonant=0,
        numeric_resonant=0,
        gamma_lims_start=[0,8], # providing starting value for eigenvalue search
        omega_lims_start=[0,3], # providing starting value for eigenvalue search
        gamma_step_start=0.2,
        omega_step_start=0.2,
        grid_factor=0.1, # smaller number → finer grid in subsequent search box
        eigenvector=false,
        progress_bar=true,
        quiet_output=false,
        quiet_scan=false,
        save=false,
        savefolder="./output/")

    if analytic_resonant==false && numeric_resonant==false
        model="non_resonant"
    elseif analytic_resonant==true
        model="analytic_resonant"
    elseif numeric_resonant==true
        model="numeric_resonant"
    end

    global nz0=gridpoints_per_pol*npol

    if save==true
        # Creating directory for storing output data
        kyscan_directory = string(savefolder,model,"_beta_",beta,"_kyrho_gene_scan_",first(kyrho_gene),"_",last(kyrho_gene),"_n_ev_",n_ev,"/")
        mkdir(betascan_directory)
    end
    
    # Creating parameters.txt file for reproducibility
    global parameters=String[]
    push!(parameters, "# KEY parameters file")
    push!(parameters, "")
    push!(parameters, "# Geometry")
    push!(parameters, "geo_name=salpha")
    push!(parameters, string("import_geo=",import_geo))
    push!(parameters, string("shear=",shear))
    push!(parameters, string("npol=",npol))
    push!(parameters, string("gridpoints=",nz0))
    push!(parameters, "")
    push!(parameters, "# Plasma")
    push!(parameters, string("beta=",beta))
    push!(parameters, string("kyrho_gene=",kyrho_gene))
    push!(parameters, string("omti=",omti))
    push!(parameters, string("omte=",omte))
    push!(parameters, string("omn=",omn))
    push!(parameters, "")
    push!(parameters, "# Simulation")
    push!(parameters, string("n_ev=",n_ev))
    push!(parameters, string("analytic_resonant=",analytic_resonant))
    push!(parameters, string("numeric_resonant=",numeric_resonant))
    push!(parameters, string("eigenvector=",eigenvector))
    push!(parameters, string("gamma_lims_start=",gamma_lims_start))
    push!(parameters, string("omega_lims_start=",omega_lims_start))
    push!(parameters, string("gamma_step_start=",gamma_step_start))
    push!(parameters, string("omega_step_start=",omega_step_start))
    push!(parameters, string("grid_factor=",grid_factor))
    push!(parameters, "")

    if save==true
        writedlm(string(betascan_directory,"parameters.txt"), parameters)
    end
 
    # Initialising matrices for saving output
    global eigenvalue_save_matrix=Complex.(zeros(length(kyrho_gene), n_ev))
    global ky_save_matrix=zeros(length(kyrho_gene), n_ev)

    if eigenvector==true
        # Initialising eigenvector-saving function 
        reset_single_M_op_efn_save() 
    end
    
    # Quiet and loud options: with and without progress bar
    progress_bar_array = ProgressBar(1:length(kyrho_gene)) 
    for_loop_array = 1:length(kyrho_gene) 

    if quiet_scan==false
        scan_array = progress_bar_array
    else
        scan_array = for_loop_array
    end

    for i in scan_array

        automatic_eigenvalue_finder_salpha(;
                                           beta=beta,
                                           shear=shear,
                                           npol=npol,
                                           gridpoints_per_pol=gridpoints_per_pol,
                                           omti=omti,
                                           omte=omte,
                                           omn=omn,
                                           import_geo=import_geo,
                                           kyrho_gene=kyrho_gene[i],
                                           n_ev=1, # number of eigenvalues to save
                                           analytic_resonant=analytic_resonant,
                                           numeric_resonant=numeric_resonant,
                                           gamma_lims_start=gamma_lims_start, 
                                           omega_lims_start=omega_lims_start, 
                                           gamma_step_start=gamma_step_start,
                                           omega_step_start=omega_step_start,
                                           grid_factor=grid_factor,
                                           progress_bar=progress_bar,
                                           quiet_output=quiet_output)
                                           
        eigenvalue_save_matrix[i,1:n_ev].=last(eigenvalue_save_converge)
            
        ky_save_matrix[i,1:n_ev].=kyrho_gene[i]
            
        if eigenvector==true
            M_op_efn_array_plot(plot=false)
            single_M_op_efn_plot(min_λ_index, plot=0) # saving eigenvector of dominant eigenmode 
        end

    end
    
    # Creating dataframe for dominant eigenvalues
    scanlog_matrix=hcat(ky_save_matrix[:,1],imag(eigenvalue_save_matrix[:,1]),real(eigenvalue_save_matrix[:,1]))
    scanlog_column_names=["kyrho_gene","gamma","omega"]
    global scanlog_dataframe=DataFrame(scanlog_matrix, scanlog_column_names)
      
    if eigenvector==true

        global eigenvector_column_names= String[]
        map(x -> push!(eigenvector_column_names, string("eigenvector",x)), 1:1:length(beta))

        # Creating dataframe for dominant eigenvectors
        eigenvector_matrix=single_M_op_efn_save_full_evecs
        global eigenvector_dataframe=DataFrame(eigenvector_matrix, eigenvector_column_names)

    end

    if save==true
        CSV.write(string(kyscan_directory,"eigenvalues.dat"),
                  scanlog_dataframe,
                  delim='\t')

        if eigenvector==true
            CSV.write(string(kyscan_directory,"eigenvectors.dat"),
                      eigenvector_dataframe,
                      delim='\t')
        end 

    end

end

##################################### End ######################################
