############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Importing GIST flux-tube geometry

function import_gist(
        fileloc::String; 
        skip=13, 
        saveplot=false, 
        gyy_rescale=false, 
        gyy_rescale_factor=1.0, 
        Ky_rescale=false, 
        Ky_rescale_factor=1.0)

    full_file = readdlm(fileloc)
    global head = full_file[1:12,:] 
    global B_ref = head[5,3]
    global my_dpdx = Float64(head[6,3])
    global q = Float64(head[7,3])
    global shat = Float64(head[8,3])
    global nz0 = Int64(head[9,3])
    global n_pol = Int64(head[10,3])
    global minor = Float64(head[4,6])
    global major = Float64(head[4,5])
    global L_ref = minor
    global ϵ_t = minor/major # inverse aspect ratio
 
    # defining θ-related quantities
    global th = range(-pi*n_pol, pi*n_pol, nz0+1)
    global th_mid = Int64(nz0/2 +1)
    global thp = th/pi

    deriv_mat_FCB() 

    # defining gist dataframe
    gist_orig = readdlm(fileloc, skipstart = skip) # obtaining matrix of gist data; skipstart removes header and 1st data point
    gist = vcat(gist_orig, transpose(gist_orig[1,:] .* [1,-1,1,1,1,1,-1,-1,1])) # appending first data point of each gist quantity to its end; ensures full nπ domain is filled; vector of ±1 ensures parity is correct (+1 for even, -1 for odd) 
    geo_names = ["gxx","gxy","gyy","B","J","Ly","Lx","dBdz","null"] # column names in dataframe
    global geo=DataFrame(gist, geo_names)

    geo[!,"Ky"] = geo.Ly .- (my_dpdx ./ (2*geo.B))
    geo[!,"shat_loc"] = dfdx_matrix_FCB * (geo.gxy ./ geo.gxx) 

    if gyy_rescale==true
        global geo.gyy = gyy_rescale_factor * geo.gyy
    end

    if Ky_rescale==true
        global geo.Ky = Ky_rescale_factor * geo.Ky
    end

    # LHS operator
    global L_op = (geo.J .* geo.B.^2).^-1 .* (dfdx_matrix_FCB * (geo.gyy ./ (geo.J .* geo.B.^2)) .* dfdx_matrix_FCB + (geo.gyy ./ (geo.J .* geo.B.^2)) .* d2fdx2_matrix_FCB) 

    print("\nYou have now imported a GIST flux-tube geometry.\n")

    geo_plot(saveplot=saveplot)

end

# Importing s-alpha analytic flux-tube geometry

function import_salpha(; 
        r = 0.18, 
        R=1, 
        q0 = 1.4, 
        shear = 0.796, 
        beta=0.0, 
        omti=1, 
        omte=0, 
        omn=1, 
        npol=5, 
        gridpoints = 5*64, 
        L_connection=false, 
        gene_implement=true, 
        saveplot=false)

    global minor = r
    global major = R
    global L_ref = major
    global B_ref = 1
    global q = q0
    global shat = shear
    global gradsum = omti + omte + (2*omn)
    global my_dpdx = beta * gradsum 
    global alpha = (q^2 * major/L_ref) * my_dpdx
    global nz0 = gridpoints
    global n_pol = npol
    global ϵ_t = minor/major

    if L_connection == true
        global L_ref = q*major
    end
    
    # defining θ-related quantities
    global th = range(-pi*n_pol, pi*n_pol, nz0+1)
    global th_mid = Int64(nz0/2 +1)
    global thp = th/pi

    deriv_mat_FCB()

    gxx = Float64.(fill(1, length(th))) 
    gxy = shat .* th .- alpha .* sin.(th)
    gyx = gxy
    gyy = 1 .+ (shat .* th .- alpha .* sin.(th)).^2 
    gyz = Float64.(fill(1/minor, length(th)))
    B = (1 .+ ϵ_t .* cos.(th)).^-1 

    if gene_implement==true
        J = abs(q/major) * (1 .+ ϵ_t .* cos.(th))
        Ly = - (B.^2) .* (L_ref/major) .* (cos.(th) .+ sin.(th) .* (shat .* th .- alpha .* sin.(th))) 
    else
        J = (1 .+ ϵ_t .* cos.(th))
        Ly = - (L_ref/major) .* (cos.(th) .+ sin.(th) .* (shat .* th .- alpha .* sin.(th))) 
    end

    Lx = - (L_ref/major) .* sin.(th) 
    dBdz = ((1 .+ ϵ_t .* cos.(th)).^-1).^2 .* ϵ_t .* sin.(th) 

    geo_mat = [gxx gxy gyy B J Ly Lx dBdz]
    geo_names = ["gxx","gxy","gyy","B","J","Ly","Lx","dBdz"] # column names in dataframe

    global geo=DataFrame(geo_mat, geo_names)

    geo[!,"Ky"] = geo.Ly .- (my_dpdx ./ (2*geo.B))
    geo[!,"shat_loc"] = dfdx_matrix_FCB * (geo.gxy ./ geo.gxx)

    # LHS operator
    global L_op = (geo.J .* geo.B.^2).^-1 .* (dfdx_matrix_FCB * (geo.gyy ./ (geo.J .* geo.B.^2)) .* dfdx_matrix_FCB + (geo.gyy ./ (geo.J .* geo.B.^2)) .* d2fdx2_matrix_FCB) 
    
    print("\nYou have now imported an analytical s-alpha flux-tube geometry.\n")

    geo_plot(saveplot=saveplot)

end

# importing s_alpha geometry from GENE output

function import_salpha_gene(fileloc::String; skip=19, saveplot=false)

    full_file= readdlm(fileloc)
    global head = full_file[1:18,:] 
    global B_ref = head[11,3]
    global my_dpdx = Float64(head[12,3])
    global q = Float64(head[3,3])
    global shat = Float64(head[4,3])
    global nz0 = Int64(head[2,3])
    global n_pol = 1
    global ϵ_t = head[8,3] # inverse aspect ratio
    global major = 1
    global minor = major*ϵ_t
    global alpha = (q^2/major)*my_dpdx 
    global L_ref = major
 
    # defining θ-related quantities
    global th = range(-pi*n_pol, pi*n_pol, nz0+1)
    global th_mid = Int64(nz0/2 +1)
    global thp = th/pi

    deriv_mat_FCB() 

    #NOTE: Jacobian here is q0/(B*major)

    # defining gist dataframe
    gist_orig = readdlm(fileloc, skipstart = skip) # obtaining matrix of gist data; skipstart removes header and 1st data point
    gist = vcat(gist_orig, transpose(gist_orig[1,:] .* [1,-1,-1,1,1,1,1,1,1,-1,1,1,1,-1,1,-1])) # appending first data point of each gist quantity to its end; ensures full nπ domain is filled; vector of ±1 ensures parity is correct (+1 for even, -1 for odd) 
    
    geo_names = ["gxx","gxy","gxz","gyy","gyz","gzz","B","dBdx","dBdy","dBdz","J","R","P","Z","dxdR","dxdZ"] # column names in dataframe
    global geo=DataFrame(gist, geo_names)

    γ1 = (geo.gxx .* geo.gyy) .- (geo.gxy .* geo.gxy) 
    γ2 = (geo.gxx .* geo.gyz) .- (geo.gxy .* geo.gxz) 
    γ3 = (geo.gxy .* geo.gyz) .- (geo.gyy .* geo.gxz) 

    geo[!,"Ly"] = (geo.dBdx .- (γ3./γ1).*geo.dBdz) # need to multiply by 1/B^2 to get same as s-alpha analytical 
    geo[!,"Ky"] = geo.Ly .- (my_dpdx ./ (2*geo.B))
    geo[!,"shat_loc"] = dfdx_matrix_FCB * (geo.gxy ./ geo.gxx)
    
    # LHS operator
    global L_op = (geo.J .* geo.B.^2).^-1 .* (dfdx_matrix_FCB * (geo.gyy ./ (geo.J .* geo.B.^2)) .* dfdx_matrix_FCB + (geo.gyy ./ (geo.J .* geo.B.^2)) .* d2fdx2_matrix_FCB) 

    print("\nYou have now imported an s-alpha geometry from GENE ouptut.\n")

    geo_plot(saveplot=saveplot)

end

function import_miller_gene(fileloc::String; skip=21, saveplot=false, gyy_rescale=false, gyy_rescale_factor=1.0, Ky_rescale=false, Ky_rescale_factor=1.0,)

    full_file = readdlm(fileloc)
    global head = full_file[1:19,:] 
    global nz0 = Int64(head[2,3])
    global q = Float64(head[3,3])
    global shat = Float64(head[4,3])
    global L_ref = head[11,3]
    global minor = Float64(head[6,2])*L_ref
    global major = Float64(head[7,2])*L_ref 
    global ϵ_t = head[8,3] # inverse aspect ratio
    global B_ref = head[12,3]
    global n_pol = head[13,3]
    global my_dpdx = Float64(head[14,3])
 
    # defining θ-related quantities
    global th = range(-pi*n_pol, pi*n_pol, nz0+1)
    global th_mid = Int64(nz0/2 +1)
    global thp = th/pi

    deriv_mat_FCB() 
    
    # defining gist dataframe
    gist_orig = readdlm(fileloc, skipstart = skip) # obtaining matrix of gist data; skipstart removes header and 1st data point
    gist = vcat(gist_orig, transpose(gist_orig[1,:] .* [1,-1,-1,1,1,1,1,1,1,-1,1,1,1,-1,1,-1])) # appending first data point of each gist quantity to its end; ensures full nπ domain is filled; vector of ±1 ensures parity is correct (+1 for even, -1 for odd) 
    
    geo_names = ["gxx","gxy","gxz","gyy","gyz","gzz","B","dBdx","dBdy","dBdz","J","R","P","Z","dxdR","dxdZ"] # column names in dataframe
    global geo=DataFrame(gist, geo_names)

    γ1 = (geo.gxx .* geo.gyy) .- (geo.gxy .* geo.gxy) 
    γ2 = (geo.gxx .* geo.gyz) .- (geo.gxy .* geo.gxz) 
    γ3 = (geo.gxy .* geo.gyz) .- (geo.gyy .* geo.gxz) 

    geo[!,"Ly"] = (geo.dBdx .- (γ3./γ1).*geo.dBdz) 
    geo[!,"Ky"] = geo.Ly .- (my_dpdx ./ (2*geo.B))
    geo[!,"shat_loc"] = dfdx_matrix_FCB * (geo.gxy ./ geo.gxx)

    if gyy_rescale==true
        global geo.gyy = gyy_rescale_factor * geo.gyy
    end

    if Ky_rescale==true
        global geo.Ky = Ky_rescale_factor * geo.Ky
    end
    
    # LHS operator
    global L_op = (geo.J .* geo.B.^2).^-1 .* (dfdx_matrix_FCB * (geo.gyy ./ (geo.J .* geo.B.^2)) .* dfdx_matrix_FCB + (geo.gyy ./ (geo.J .* geo.B.^2)) .* d2fdx2_matrix_FCB) 

    print("\nYou have now imported a miller geometry from GENE ouptut.\n")

    geo_plot(saveplot=saveplot)

end

##################################### End ######################################
