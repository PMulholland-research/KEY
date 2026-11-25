# KEY — Kinetic-Ballooning-Mode Eigenvalue Yielder

**KEY** (Kinetic-Ballooning-Mode Eigenvalue Yielder) is a Julia codebase that extracts linear growth rates, frequencies, and eigenfunctions from gyrokinetic- or fluid-type equations. It provides a partial description of plasma performance in stellarator and tokamak magnetic-confinement fusion reactors. 

Taking flux-tube geometry and kinetic profile information as input, KEY uses an iterative refinement process to yield eigenmode information for ion-temperature-gradient-driven (ITG) and kinetic-ballooning-type (KBM) microinstabilities.

> **Reference:** The underlying theory and first published results of KEY can be found in:  
> *Journal of Plasma Physics*, Volume 91, Issue 5, October 2025, E142  
> DOI: [10.1017/S002237782510086X](https://doi.org/10.1017/S002237782510086X)

---

## Table of Contents
- [Prerequisites](#prerequisites)
- [Getting Started](#getting-started)
- [Usage](#usage)
  - [Importing Flux-Tube Geometry](#importing-flux-tube-geometry)
  - [Calculating Eigenvalues and Eigenvectors](#calculating-eigenvalues-and-eigenvectors)
  - [Parameter Scans](#parameter-scans)
- [License & Copyright](#license--copyright)

---

## Prerequisites

To run KEY, install the following required Julia packages:

- `CSV`
- `DataFrames`
- `DelimitedFiles`
- `FastGaussQuadrature`
- `GLMakie`
- `Interpolations`
- `KissSmoothing`
- `LaTeXStrings`
- `LinearAlgebra`
- `ProgressBars`
- `SpecialFunctions`
- `Trapz`

<details>
<summary>Quick install snippet (click to expand)</summary>

```julia
using Pkg
Pkg.add([
    "CSV", "DataFrames", "DelimitedFiles", "FastGaussQuadrature",
    "GLMakie", "Interpolations", "KissSmoothing", "LaTeXStrings",
    "LinearAlgebra", "ProgressBars", "SpecialFunctions", "Trapz"
])
```
</details>

---

## Getting Started

1. Open a terminal and launch Julia in the directory where `main_KEY.jl` is located.
2. Include the main script to load the codebase:

```julia
include("main_KEY.jl")
```

---

## Usage

### Importing Flux-Tube Geometry

* **Import a GIST geometry file:**
  ```julia
  import_gist("path/to/gist_file.dat")
  ```
* **Import an $s$-$\alpha$ tokamak geometry:**
  ```julia
  import_salpha()
  ```

> **Note:** For more information on these functions, see [`src/import_geo.jl`](src/import_geo.jl).

---

### Calculating Eigenvalues and Eigenvectors

* **Evaluate the eigenvalue plane $(\gamma, \omega_r)$** for a fixed set of input parameters:
  ```julia
  condition_M_op(beta)
  ```
* **Calculate the eigenfunction** for a given eigenvalue input:
  ```julia
  M_op_efn_array_plot_any_input(eigenvalue)
  ```

> **Note:** For more information on these functions, see [`src/matrix_formalism.jl`](src/matrix_formalism.jl).

---

### Parameter Scans

* **Perform an automatic $\beta$-scan** using GIST geometry:
  ```julia
  automatic_betascan_general_geometry_zoomer()
  ```
* **Perform an automatic $k_y$-scan** using GIST geometry:
  ```julia
  automatic_kyscan_general_geometry_zoomer()
  ```

> **Notes:** 
> - Analogous functions exist for $s$-$\alpha$ tokamak geometry, whereby the geometry is generated consistently throughout the scan.
> - An input-parameters list and eigenvalue-eigenvector data can be automatically stored to the output directory using these functions.
> - For full details, see [`src/automatic_parameter_scan.jl`](src/automatic_parameter_scan.jl).

---

## License & Copyright

**KEY — Kinetic-Ballooning-Mode Eigenvalue Yielder**  
Copyright © 2026 Stichting Nederlandse Wetenschappelijk Onderzoek Instituten (NWO-I), acting through its institute DIFFER (Dutch Institute for Fundamental Energy Research), and Technische Universiteit Eindhoven. Both parties hold the copyright jointly in equal shares.

**Original Author:** Paul Mulholland

KEY is free software: you can redistribute it and/or modify it under the terms of the **GNU Lesser General Public License** as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

KEY is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License and the GNU General Public License along with KEY (files [`COPYING.LESSER`](COPYING.LESSER) and [`COPYING`](COPYING)). If not, see <https://www.gnu.org/licenses/>.