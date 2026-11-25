############################ Software license info #############################

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.

#################################### Begin #####################################

# Required Julia packages

using CSV
using DataFrames 
using DelimitedFiles 
using FastGaussQuadrature 
using GLMakie 
using Interpolations 
using KissSmoothing 
using LaTeXStrings 
using LinearAlgebra 
using ProgressBars 
using SpecialFunctions 
using Trapz 

# Importing source-code scripts

include("./src/automatic_parameter_scan.jl")
include("./src/finite_diff_stencil.jl")
include("./src/import_geo.jl")
include("./src/kbm_equation_RHS.jl")
include("./src/matrix_formalism.jl")
include("./src/misc_functions.jl")
include("./src/plotting_functions.jl")

# Introductory text for the user

print("\n############################ Software license info #############################
\n# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2026 NWO-I (DIFFER) and Technische Universiteit Eindhoven
#
# This file is part of KEY (Kinetic-Ballooning-Mode Eigenvalue Yielder).
# KEY is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License, version 3 or (at your
# option) any later version. See the files COPYING.LESSER and COPYING.\n")

print("\n############### Welcome to the KBM-eigenvalue-yielder (KEY) code ###############
\n# KEY generates KBM eigenvalues and eigenvectors based on flux-tube geometry and 
# plasma parameters.
#
# To begin, you can do one of the following:
# i) Import the in-built analytic s-alpha (circular tokamak) geometry.
# ii) Import a flux-tube geometry file (e.g., GIST).
#
# Following this, you are ready to compute KBM eigenvalues and eigenvectors.
\n################################################################################\n")

##################################### End ######################################
