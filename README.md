# Boundary self-energy chain closure

MATLAB code accompanying the manuscript **Boundary Self-Energy Calculus for Chain Representations of Open Quantum Systems**, by Pengju Chen, Da-wei Luo, and Ting Yu.

This package compares boundary self-energy (BSE) closure, hard termination, and a complex absorbing potential (CAP) for a chain representation of an Ohmic environment coupled to a qubit. It computes the qubit excited-state population and the broadened effective spectrum at the qubit port. The default benchmark uses a **300-site chain as the reference in both the time and frequency domains**.

The manuscript establishes the equivalence of the star-geometry and exact chain representations. The numerical comparisons in this package retain the finite-chain reference.

## Requirements

Use MATLAB with `exportgraphics`, which was introduced in R2020a. No third-party packages or specialized toolbox calls were identified in the supplied source. Compatibility with a particular MATLAB release must be established by running the checks below; GNU Octave compatibility has not been tested.

See the [MathWorks documentation for exportgraphics](https://www.mathworks.com/help/matlab/ref/exportgraphics.html).

## Quick start

Place all 12 MATLAB files in the same folder. Open that folder as the MATLAB current folder, then run:

```matlab
verification = verify_bse_solver();
writetable(verification,'solver_verification.csv');
run_chain_closure_benchmark_fixed
```

The benchmark generates five numbered figures, one diagnostic figure, and three summary CSV files in the current folder. Each figure is exported as a vector PDF and a 300 dpi PNG. Running the script again overwrites outputs with the same names.

**To archive the full numerical curves**, execute the save command in [REPRODUCIBILITY.md](REPRODUCIBILITY.md) immediately after the benchmark. The benchmark itself exports figures and summary tables, but does not automatically save the full time and frequency arrays.

The main benchmark and auxiliary scan scripts clear the workspace and close existing figures. The main benchmark also sets session graphics defaults. Run each script as a complete script from this package folder.

## Default comparison

The default parameters are `omega_c = 1`, `alpha = 0.02`, and `epsilon = 0.5`, with `t = 0:0.1:300`. The retained chain lengths are `N = [10, 25, 40, 60]`. All CAP comparisons use the same coefficients, `aCAP = 0` and `bCAP = 0.57`, with onset indices `[5, 12, 20, 30]`, respectively.

The main time-domain and spectrum figures compare the reference with hard termination at `N = 60`, CAP at `(N, Nd) = (10, 5)` and `(60, 30)`, and BSE closure at `N = 10`. An additional BSE scan over `N = 1:60` determines the operational closure depth for three error tolerances.

The reference-length rule gives `Nref = 300` for the supplied defaults. It can increase the reference length if the observation window is extended. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the precise rule and the BSE indexing convention.

## Files

| File | Purpose |
| --- | --- |
| `run_chain_closure_benchmark_fixed.m` | Main benchmark, figure generation, and summary CSV export |
| `scan_cap_coefficients.m` | Coarse and fine CAP parameter scans using the population-error objective |
| `verify_bse_solver.m` | Solver checks against a finite chain with the same BSE closure model |
| `check_bse_kgrid_convergence.m` | Separate long-time quadrature convergence check |
| `woods_chain_coeffs.m` | Ohmic chain coefficients and qubit coupling |
| `boundary_self_energy_chain.m` | Inward BSE recurrence for a specified chain and boundary condition |
| `cap_profile.m` | Linear plus quadratic absorption profile |
| `finite_chain_amplitude_spectral.m` | Finite-chain amplitude from eigendecomposition |
| `sigma_uniform_tail.m` | Homogeneous semi-infinite tail self-energy and its derivative |
| `green_bse.m` | Qubit Green function for the matched BSE closure |
| `find_bound_states_from_green.m` | Bound-state energies and residues outside the continuum band |
| `bse_amplitude_kquad.m` | Continuum quadrature and bound-state contributions to the amplitude |
| [REPRODUCIBILITY.md](REPRODUCIBILITY.md) | Parameters, methods, verification commands, and data archiving |
| [DATA_DICTIONARY.md](DATA_DICTIONARY.md) | Output filenames, table columns, and saved-array meanings |
| [CITATION.bib](CITATION.bib) | Manuscript citation |

## Additional runs

```matlab
scan_cap_coefficients
```

This reruns CAP selection and writes two CSV files. It does not automatically modify the CAP coefficients in the main benchmark.

```matlab
check_bse_kgrid_convergence
```

This runs a separate long-time convergence check with different physical parameters and displays a diagnostic plot. It is not the default manuscript benchmark.

## Data and citation

The scripts generate the numerical results locally. Consult [DATA_DICTIONARY.md](DATA_DICTIONARY.md) for the exported data and the optional full-curve archive. When releasing the repository, include the generated data corresponding to the manuscript, along with the MATLAB version and verification output.

If you use this work, cite the accompanying manuscript using [CITATION.bib](CITATION.bib). The entry currently identifies an unpublished manuscript; update it when publication metadata become available.

For scientific questions, contact Pengju Chen at pchen9@stevens.edu. For reproducibility issues, include the MATLAB version, the script used, any parameter changes, and the complete error message in a repository issue.

## License

A software license has not yet been specified for this package.
