############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Definitions
#
# Magnetic field strength:
# -> B(θ): field-line dependent magnetic-field strength
# -> B_ref: reference magnetic-field strength on axis
# -> B^hat ≡ B(θ)/B_ref
#
# Normalized plasma pressure:
# -> β_ref ≡ 8πi * n_e0 * T_e0 / B_ref^2 
# -> β_e ≡ 8πi n_e0 T_e0 / B(θ)^2 
# -> β_i ≡ 8πi n_i0 T_i0 / B(θ)^2 
# -> β_e = β_ref / (B^hat)^2 
#
# For n_e0 ≈ n_i0, and τ ≡ T_e0/T_i0 we have
# -> β_e = τ * β_i
#
# In KEY, β = β_ref (as used in GENE), which depends on B_ref^-2 and not on B(θ)^-2
#
# Important distinction between KEY and GENE: 
# -> ρ_ref^GENE ≡ ρ_s
# -> ρ_ref^KEY = √2 * ρ_ref^GENE = √2 * ρ_s 
#
# This means that kyrho used here in KEY is NOT the same as kyrho in GENE 
#
# Therefore, have created a converter below for k_y*ρ_ref^KEY and k_y*ρ_ref^GENE, 
# as the user expects to insert kyrho as defined in GENE
#
# Throughout this script, the argument 'kyrho' refers to kyrho_KEY 

##################### Defining quantities in KBM equation ######################

# Converting from kyrho^GENE to kyrho_KEY
kyrho_key(kyrho_gene) = √2 * kyrho_gene 

# Density-gradient scale length
L_n(omn) = L_ref/omn # for L_ref = a (R), omn = a/L_n (R/L_n)

# Normalized inverse density gradient
ϵ_n(omn) = L_n(omn)/L_ref 

# α_MHD (available after importing s-alpha geometry)
alpha_MHD(β, omti, omte, omn) = (q^2 * major/L_ref) * β * (omti + omte + 2*omn)  

# Alfvén-frequency prefactor (analogous to term in Cheng, PoF, 1982) 
# Ω_A^-2 = (β_i * τ) / (2 * ϵ_n^2) 
#        = β_e / (2 * ϵ_n^2) 
#        = β_ref / (2 * (B^hat)^2 * ϵ_n^2)
Ω_A(β,τ,omn) = (β ./ (2 .* geo.B.^2 .* ϵ_n(omn)^2)).^-0.5 

# FLR term
# b = 0.5 * (k_⟂ * ρ_thi)^2 
#   ≈ 0.5 * (k_y * ρ_thi)^2 * g^yy (k_x = 0)
#   = 0.5 * (k_y * ρ_ref/B^hat)^2 * g^yy 
b(kyrho) = 0.5 .* (kyrho ./ geo.B).^2 .* geo.gyy 

# Curvature terms

# Defining grad-B and curvature drift frequencies  

Ω_B(omn,τ) = (ϵ_n(omn)/τ) * (geo.Ly ./ geo.B) # Ω_B = ω_B/ω_*e 
Ω_curv(omn,τ) = (ϵ_n(omn)/τ) * (geo.Ky ./ geo.B) # Ω_κ = ω_κ/ω_*e 

# Defining different options for magnetic-curvature-drift term (pressure term in GENE): Full-drift, gradB_eq_curv, curv_eq_gradB

# Full drift: ω_κ = ω_B - dpdx/2B
Ω_fulldrift(omn,τ) = Ω_B(omn,τ) + Ω_curv(omn,τ)

# gradB_eq_curv drift: ω_B → ω_κ such that (ω_B + ω_κ) ≡ 2ω_κ
Ω_gradBeqcurv(omn,τ) = 2*Ω_curv(omn,τ) # grad-B drift = curvature drift 

# curv_eq_gradB drift: ω_κ → ω_B such that (ω_B + ω_κ) ≡ 2ω_B
Ω_curveqgradB(omn,τ) = 2*Ω_B(omn,τ) # curvature drift = grad-B drift 

# Choosing magnetic-curvature-drift 
# NOTE: Be careful with this; might replace this term with separate grad-B and curvature frequencies 
Ω_mag_curv(omn,τ) = Ω_gradBeqcurv(omn,τ) # (Ω_B + Ω_κ) ≡ 2Ω_κ # see notes for explanation of choice 

# Sign/direction convention in KEY: 
# -> Ion direction: negative
# -> Electron direction: positive
#
# Normalized diamagnetic drift frequency:
# -> ω_*i / (v_thi/L_ref) = - 1/2 * kyrho * L_ref/L_n
#
# Using 
# -> v_thi = sqrt(2*T_i0/m_i) 
# -> c_s = sqrt(T_e0/m_i) 
# -> τ ≡ T_e0/T_i0
#
# we have
# -> v_thi = sqrt(2/τ)*c_s: 
#
# such that the normalized diamagnetic drift frequency is:
# -> ω_star_i_norm ≡ ω_*i / (c_s/L_ref) = - sqrt(2/τ) * 1/2 * kyrho * L_ref/L_n 

# Normalized diamagnetic drift frequencies:
ω_star_i_norm(kyrho,omn,τ) = - sqrt(2/τ) * 0.5 * kyrho * omn
ω_star_e_norm(kyrho,omn,τ) = - τ * ω_star_i_norm(kyrho,omn,τ)

α_i(kyrho,omn,n,Ω,η_i,τ) = 1 + (1/(τ*Ω))*(1 + n*η_i) 
α_e(kyrho,omn,n,Ω,η_e) = 1 - (1/Ω)*(1 + n*η_e) 
ω_p_i(kyrho,omn,τ,η_i,η_e) = ω_star_i_norm(kyrho,omn,τ)*(1 + η_i)
ω_p_e(kyrho,omn,τ,η_i,η_e) = - ω_star_e_norm(kyrho,omn,τ)*(1 + η_e)
ω_p(kyrho,omn,τ,η_i,η_e) = ω_p_i(kyrho,omn,τ,η_i,η_e) + ω_p_e(kyrho,omn,τ,η_i,η_e)

Ω_p_i(τ,η_i) = -(1/τ)*(1+η_i)

# Using
# ω_norm = ω / (c_s/L_ref)
# ω_star_e_norm = ω_*e / (c_s/L_ref)

# we have 
# -> Ω ≡ ω / ω_*e = ω_norm / ω_star_e_norm 
# -> ω ≡ Ω * ω_*e ⟹  ω_norm = Ω * ω_star_e_norm

Ω(ω_norm,kyrho,omn,τ) = ω_norm / ω_star_e_norm(kyrho,omn,τ) # converts GENE ω_norm to KEY Ω (sign is reversed: KEY assumes electron direction is positive)
ω_norm(Ω,kyrho,omn,τ) = Ω * ω_star_e_norm(kyrho,omn,τ) # converts KEY Ω to GENE ω_norm (but sign is reversed: KEY assumes electron direction is positive) 

# Gamma function Γ - needed for analytic-resonant approach
Γ_orig(n,τ,kyrho) = besseli.(n, b(kyrho)) .* exp.(-b(kyrho)) 

# to combat overflows at high b(kyrho), need to put cap on b(kyrho) at ~700
function Γ(n,τ,kyrho; threshold = 700)
    
    global gamfn_ar = Vector{Float64}(undef, nz0+1) 
    
    gamfn_ar[b(kyrho) .< threshold] = besseli.(n, b(kyrho)[b(kyrho) .< threshold]) .* exp.(-b(kyrho)[b(kyrho) .< threshold]) 
    gamfn_ar[b(kyrho) .> threshold] .= besseli(n, 700) * exp(-700) 

    return gamfn_ar

end

################################ Ion integral Q ################################

# Defined for different frameworks:
# i)    Q_fluid (non-resonant)
# ii)   Q_analytic_resonant (includes partial toroidal resonance; assumes v_par^2 ≈ 0.5*v_perp^2, gradB_eq_curv drift) 
# iii)  Q_numeric_resonant (includes full toroidal resonance; can be used with full-drift, gradB_eq_curv drift or curv_eq_gradB drift))

# i) Q_fluid
Q_fluid(Ω,τ,kyrho,η_i,omn) = α_i(kyrho,omn,0,Ω,η_i,τ) .+ ((Ω_mag_curv(omn,τ)/Ω) .- b(kyrho)) .* α_i(kyrho,omn,1,Ω,η_i,τ)   

# ii) Q_kinetic_partial_res

# Using trick by Maikel Morren, where we transform:
# ξ → ξ_new = sign(imag[ξ]) * ξ
# which ensures that imag[ξ_new] > 0, i.e., that Q_cheng is valid everywhere.

# In short, we know that Q_integral has complex-conjugate symmetry built in:
# → Q^*(Ω^*) ≡ Q(Ω)

# However, Cheng rewrites this Q_integral using the plasma dispersion function Z(ξ),
# which has limited validity, i.e., it only return valid results when imag[ξ] > 0,
# which only occurs when sign[γ] = sign[κ], i.e., when growth rates and curvature have the same sign.
 
# Maikel's transform restores the complex-conjugate symmetry to Q_cheng,
# while still retaining the nice analytic properties of Z(ξ).

# Plasma dispersion function Z
ξ(Ω,τ,omn) = sqrt.(Complex.( Ω ./ (2*Ω_mag_curv(omn,τ)) )) 
Z(Ω,τ,omn) = im*sqrt(π) .* exp.(-ξ(Ω,τ,omn).^2) .* (1 .+ im*erfi.(sign.(imag(ξ(Ω,τ,omn))).*ξ(Ω,τ,omn))) 

# D_0 and D_1
D_0(Ω,τ,omn) = ((-sign.(imag(ξ(Ω,τ,omn))).*ξ(Ω,τ,omn))/Ω) .* Z(Ω,τ,omn) 
D_1(Ω,τ,omn) = ((-ξ(Ω,τ,omn).^2)/Ω) .* (1 .+ sign.(imag(ξ(Ω,τ,omn))).*ξ(Ω,τ,omn) .* Z(Ω,τ,omn)) 

# Resonant ion integral Q
# Note: τ^1 in front because Cheng uses Resonant-integral ≡ τQ, not Q
Q_cheng(Ω,τ,kyrho,η_i,omn) = (τ^-1) .* ( ((τ*Ω + 1 - (η_i/2)) .* Γ(0,τ,kyrho) .+ η_i .* b(kyrho) .* (Γ(1,τ,kyrho) .- Γ(0,τ,kyrho))) .* D_0(Ω,τ,omn) .+ (η_i .* Γ(0,τ,kyrho) .* D_1(Ω,τ,omn)) ) 

# Alternative expression for Q_cheng which is valid for Ω_mag_curv→0
Q_nocurv(Ω,τ,kyrho,η_i,omn) = (τ*Ω)^-1 .* ( (η_i .* Γ(0,τ,kyrho))./2 .+ (τ*Ω + 1 - η_i/2) .* Γ(0,τ,kyrho) .+ (η_i .* b(kyrho)) .* (Γ(1,τ,kyrho) .- Γ(0,τ,kyrho)))  

# Q_cheng breaks down when Ω_mag_curv→0;
# The following function creates a hybrid-Q array, composed of Q_cheng (where applicable) and Q_nocurv/Q_fluid, 
# where the latter are appropriate when Ω_mag_curv→0 or when Ω_mag_curv << Ω

function Q_analytic_resonant(Ω,τ,kyrho,η_i,omn; threshold = 5) # threshold is hand-tuned 

    Q_ar = Vector{ComplexF64}(undef, nz0+1)

    Q_ar[abs.(real(ξ(Ω,τ,omn))) .< threshold .&& abs.(imag(ξ(Ω,τ,omn))) .< threshold] = Q_cheng(Ω,τ,kyrho,η_i,omn)[abs.(real(ξ(Ω,τ,omn))) .< threshold .&& abs.(imag(ξ(Ω,τ,omn))) .< threshold]
    Q_ar[abs.(real(ξ(Ω,τ,omn))) .> threshold .|| abs.(imag(ξ(Ω,τ,omn))) .> threshold] = Q_fluid(Ω,τ,kyrho,η_i,omn)[abs.(real(ξ(Ω,τ,omn))) .> threshold .|| abs.(imag(ξ(Ω,τ,omn))) .> threshold]

    return Q_ar

end

# iii)  Q_numeric_resonant
# Numerical evaluation of resonant ion integral Q
# See notes for explanation of steps
    
# Defining different options for magnetic-curvature-drift term (pressure term in GENE): Full-drift, gradB_eq_curv, curv_eq_gradB

# Full toroidal resonance integral

Q_integrand_full_toroidal_resonance_full_drift(Ω,τ,kyrho,η_i,omn,μ,ν) = besselj0.(sqrt.(2*b(kyrho)*μ)).^2 .* (τ*Ω + 1 + η_i*(ν^2 + μ - 3/2)) .* (Ω .- 2*Ω_curv(omn,τ)*ν^2 .- Ω_B(omn,τ)*μ).^-1

Q_integrand_full_toroidal_resonance_gradBeqcurv_drift(Ω,τ,kyrho,η_i,omn,μ,ν) = besselj0.(sqrt.(2*b(kyrho)*μ)).^2 .* (τ*Ω + 1 + η_i*(ν^2 + μ - 3/2)) .* (Ω .- 2*Ω_curv(omn,τ)*ν^2 .- Ω_curv(omn,τ)*μ).^-1

# Curvature model (Cheng) 

Q_integrand_curvature_model_gradBeqcurv_drift(Ω,τ,kyrho,η_i,omn,μ,ν) = besselj0.(sqrt.(2*b(kyrho)*μ)).^2 .* (τ*Ω + 1 + η_i*(ν^2 + μ - 3/2)) .* (Ω .- 4*Ω_curv(omn,τ)*ν^2).^-1

Q_integrand_curvature_model_curveqgradB_drift(Ω,τ,kyrho,η_i,omn,μ,ν) = besselj0.(sqrt.(2*b(kyrho)*μ)).^2 .* (τ*Ω + 1 + η_i*(ν^2 + μ - 3/2)) .* (Ω .- 4*Ω_B(omn,τ)*ν^2).^-1

# grad-B model (Terry) 

Q_integrand_gradB_model_gradBeqcurv_drift(Ω,τ,kyrho,η_i,omn,μ,ν) = besselj0.(sqrt.(2*b(kyrho)*μ)).^2 .* (τ*Ω + 1 + η_i*(ν^2 + μ - 3/2)) .* (Ω .- 2*Ω_curv(omn,τ)*μ).^-1

Q_integrand_gradB_model_cuveqgradB_drift(Ω,τ,kyrho,η_i,omn,μ,ν) = besselj0.(sqrt.(2*b(kyrho)*μ)).^2 .* (τ*Ω + 1 + η_i*(ν^2 + μ - 3/2)) .* (Ω .- 2*Ω_B(omn,τ)*μ).^-1

# Choosing magnetic-curvature-drift options
# NOTE: for all options, remember to change Ω_mag_curv above!

#Q_integrand(Ω,τ,kyrho,η_i,omn,μ,ν) = Q_integrand_full_toroidal_resonance_full_drift(Ω,τ,kyrho,η_i,omn,μ,ν) # Don't use unless δB_∥ is implemented 
Q_integrand(Ω,τ,kyrho,η_i,omn,μ,ν) = Q_integrand_full_toroidal_resonance_gradBeqcurv_drift(Ω,τ,kyrho,η_i,omn,μ,ν) # Use this! 

#Q_integrand(Ω,τ,kyrho,η_i,omn,μ,ν) = Q_integrand_curvature_model_gradBeqcurv_drift(Ω,τ,kyrho,η_i,omn,μ,ν)  
#Q_integrand(Ω,τ,kyrho,η_i,omn,μ,ν) = Q_integrand_curvature_model_curveqgradB_drift(Ω,τ,kyrho,η_i,omn,μ,ν) 

#Q_integrand(Ω,τ,kyrho,η_i,omn,μ,ν) = Q_integrand_gradB_model_gradBeqcurv_drift(Ω,τ,kyrho,η_i,omn,μ,ν) 
#Q_integrand(Ω,τ,kyrho,η_i,omn,μ,ν) = Q_integrand_gradB_model_cuveqgradB_drift(Ω,τ,kyrho,η_i,omn,μ,ν) 

# Predefining Gaussian-quadrature nodes
# μ: v_perp^2
# ν: v_par

# Near marginality
#Nmu=21 # use 21
#Nnu=31 # use 31 

# Away from marginality
Nmu=11 # use 11
Nnu=21 # use 21 

xmu, wmu = gausslaguerre(Nmu)
xnu, wnu = gausshermite(Nnu)

weight_mat = wmu .* transpose(wnu)

mumat = repeat(xmu, 1, Nnu)
numat = repeat(transpose(xnu), Nmu, 1)

function Q_numeric_resonant(Ω,τ,kyrho,η_i,omn); 

    # Filling Q-integrand matrix
    global Q_integrand_mat = Q_integrand.(Ω,τ,kyrho,η_i,omn,mumat,numat)

    return (τ*sqrt(π))^-1 * sum(weight_mat .* Q_integrand_mat) 

end

############################# RHS of KBM equation ##############################

# RHS: fluid 
RHS_fluid(Ω,τ,kyrho,η_i,η_e,omn,β) = -Ω_A(β,τ,omn).^(-2) .* ( -(τ^(-1) * (1+η_i) + (1+η_e) ) * τ .* Ω_mag_curv(omn,τ) .+ b(kyrho) * Ω * (Ω + (1+η_i)/τ) ) 

# RHS: analytic_resonant
RHS_analytic_resonant(Ω,τ,kyrho,η_i,η_e,omn,β) = Ω^2 .* (Ω_A(β,τ,omn).^-2) .* τ .* ( (Q_analytic_resonant(Ω,τ,kyrho,η_i,omn) .- α_i(kyrho,omn,0,Ω,η_i,τ)) .* α_e(kyrho,omn,0,Ω,η_e) ./ (1 + τ .- τ*Q_analytic_resonant(Ω,τ,kyrho,η_i,omn)) .- α_e(kyrho,omn,1,Ω,η_e) .* (Ω_mag_curv(omn,τ)/Ω) ) 

# RHS: numeric_resonant
RHS_numeric_resonant(Ω,τ,kyrho,η_i,η_e,omn,β) = Ω^2 .* (Ω_A(β,τ,omn).^-2) .* τ .* ( (Q_numeric_resonant(Ω,τ,kyrho,η_i,omn) .- α_i(kyrho,omn,0,Ω,η_i,τ)) .* α_e(kyrho,omn,0,Ω,η_e) ./ (1 + τ .- τ*Q_numeric_resonant(Ω,τ,kyrho,η_i,omn)) .- α_e(kyrho,omn,1,Ω,η_e) .* (Ω_mag_curv(omn,τ)/Ω) )  

##################################### End ######################################
