# FastDysonSolver

**FastDysonSolver** is a reproducible scientific project implemented in the [Julia Language](https://julialang.org/) and structured with [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/). This repository accompanies the research paper [Solving the Transient Dyson Equation with Quasilinear Complexity via Matrix Compression](https://doi.org/10.1103/q222-w14m) by *Baptiste Lamic*
This repository contains the code required to *reproduce the figures* presented in the manuscript. The *core solver* is archived on *Zenodo* [https://zenodo.org/badge/DOI/10.5281/zenodo.17560811.svg](https://doi.org/10.5281/zenodo.17560811) and the source code is hosted on GitHub [github](https://github.com/BaptisteLamic/NonEquilibriumGreenFunction.jl/tree/0.2.5). 


## Getting Started

To reproduce the results locally:

1. **Clone this repository.**
2. **Open a Julia REPL and run:**
   ```julia
   using Pkg
   Pkg.add("DrWatson") # Required for project activation
   Pkg.activate(".")
   Pkg.instantiate()
   ```
   This will install all necessary dependencies and set up the project environment.

Most scripts begin with:
```julia
using DrWatson
@quickactivate "FastDysonSolver"
```
This ensures the project is activated and paths are handled correctly.

## Reproducing Figures

All figures from the paper can be generated using the provided Jupyter notebooks:

- [SQDS_junction_benchmarks.ipynb](notebooks/SQDS_junction_final.ipynb)
- [SQDS_junction_evaluate_steady_state.ipynb](notebooks/SQDS_junction_evaluate_steady_state.ipynb)

The original computations were performed on an Intel CPU Ultra 265k with 48 GB RAM. However, with appropriate parameter tuning, most calculations can be reproduced with 16 GB RAM.

## Data and Results

Simulation results for the IV curve (Figure 3 in the paper) are available in the [data/simulations/steady_state_current](data/simulations/steady_state_current) directory. These simulations are computationally intensive and may take several hours to run. However, all results are cached and reused if available, so you can generate the corresponding plots directly using the provided data without rerunning the simulations:

- Run [SQDS_junction_evaluate_steady_state.ipynb](notebooks/SQDS_junction_evaluate_steady_state.ipynb) to reproduce the IV curve plots.

## Extra

The notebook [SQDS_referee](notebooks/SQDS_referee.ipynb) performs simulations analogous to those presented in [Cheng et al., Phys. Rev. B 110, 125417 (2024)](https://doi.org/10.1103/PhysRevB.110.125417).
While the corresponding figure is not included in the manuscript, it was used to clarify discussions with one of the referees. 

## License

This project is licensed under the GNU General Public License v3.0.  
See the [LICENSE](LICENSE) file for details.

## Author

Baptiste Lamic

---

For questions or issues, please open an issue on this repository.