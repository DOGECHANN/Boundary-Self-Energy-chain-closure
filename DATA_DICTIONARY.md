# Output data dictionary

All filenames below refer to files in the MATLAB current folder. CSV files contain a header row and no row-name column. Results correspond to the parameters in the source revision used to generate them.

The supplied benchmark automatically exports figures and summary CSVs. Full numerical curves are saved only when the additional `save` command in [REPRODUCIBILITY.md](REPRODUCIBILITY.md) is executed.

## Main benchmark CSV files

### `termination_results.csv`

Four rows under the default settings, in the order `N = 10, 25, 40, 60`.

| Column | Meaning |
| --- | --- |
| `N` | Compared chain length or matched BSE closure depth, using the convention described in REPRODUCIBILITY.md |
| `CAPOnsetNd` | CAP onset index paired with `N` |
| `HardTimeMaxError` | Maximum sampled population error for hard termination |
| `CAPTimeMaxError` | Maximum sampled population error for CAP |
| `BSETimeMaxError` | Maximum sampled population error for BSE closure |
| `HardSpectralL2` | Relative spectral L2 error for hard termination |
| `CAPSpectralL2` | Relative spectral L2 error for CAP |
| `BSESpectralL2` | Relative spectral L2 error for BSE closure |
| `BSENk` | Number of wave-number quadrature nodes actually used |
| `BSESpectralWeight` | Continuum weight plus bound-state weights |
| `BSENormalizationError` | Absolute difference between total spectral weight and one |
| `BSEBoundWeight` | Sum of bound-state weights |

All six comparison errors use the same 300-site reference under the default settings. The spectral errors use the common broadening `etaSpec = 2.5/Nref`. Population errors, relative spectral errors, and spectral weights are dimensionless.

### `operational_closure_scan.csv`

Sixty rows under the default settings.

| Column | Meaning |
| --- | --- |
| `N` | BSE closure depth, from 1 through 60 |
| `BSETimeMaxError` | Maximum sampled population error relative to the reference |

### `operational_closure_depths.csv`

Three rows, one for each tolerance.

| Column | Meaning |
| --- | --- |
| `Tolerance` | Required upper bound on the sampled population error |
| `OperationalClosureDepth` | Smallest tested depth satisfying the tolerance |

If no tested depth qualifies, MATLAB stores `NaN`. Depending on MATLAB's CSV serialization, a missing value may appear as an empty field; interpret it as no qualifying depth within the scanned range, not as zero.

## CAP scan CSV files

These are produced by `scan_cap_coefficients.m`.

### `cap_parameter_scan.csv`

| Column | Meaning |
| --- | --- |
| `aCAP` | Linear absorption-rate coefficient |
| `bCAP` | Quadratic absorption-rate coefficient |
| `GeometricMeanTimeError` | Geometric mean of the four sampled maximum population errors |
| `TimeErrorN10` | Population error for `N = 10`, `Nd = 5` |
| `TimeErrorN25` | Population error for `N = 25`, `Nd = 12` |
| `TimeErrorN40` | Population error for `N = 40`, `Nd = 20` |
| `TimeErrorN60` | Population error for `N = 60`, `Nd = 30` |

The first 165 rows are the coarse grid; the remaining rows are the fine grid. There is no scan-stage column. Parameter pairs can appear in both stages. The error score uses `exp(mean(log(max(errors,realmin))))`.

### `cap_selected_profile_results.csv`

Four rows for the selected profile after rounding its coefficients to two decimal places.

| Column | Meaning |
| --- | --- |
| `N` | Chain length |
| `CAPOnsetNd` | CAP onset index |
| `CAPTimeMaxError` | Maximum sampled population error |
| `CAPSpectralL2` | Relative broadened-spectrum L2 error |

The selected coefficients are printed by the script. Keep the run log with this table, since the table does not contain columns for those coefficients.

## Optional verification CSV

`solver_verification.csv` is created by the explicit `writetable` command in the reproducibility instructions, not automatically by `verify_bse_solver.m`.

| Column | Meaning |
| --- | --- |
| `alpha` | Bath coupling parameter |
| `epsilon` | Qubit frequency |
| `N` | BSE closure depth |
| `TotalWeight` | Continuum plus bound-state spectral weight |
| `BoundStateCount` | Number of bound states found |
| `AmplitudeErrorVsSameClosure` | Maximum amplitude difference against a finite chain implementing the same closure |
| `PopulationErrorVsSameClosure` | Corresponding maximum population difference |

There are 12 rows for the three parameter cases and four closure depths. These errors assess the numerical solver for the specified closure, rather than its approximation to the original mapped bath.

## Optional `benchmark_data.mat`

The save command in REPRODUCIBILITY.md produces this file.

| Variable(s) | Meaning |
| --- | --- |
| `omega_c`, `alpha`, `epsilon`, `g` | Physical parameters and qubit-to-chain coupling |
| `omegaInf`, `kappaInf` | Homogeneous tail frequency and coupling |
| `omega`, `kappa` | Mapped bath frequencies and couplings; MATLAB index `n` denotes bath site/link `n` |
| `t`, `tMax`, `dt` | Time grid, end time, and spacing |
| `Nref` | Reference bath length |
| `Nscan`, `NdScan` | Compared depths and corresponding CAP onset indices |
| `aCAP`, `bCAP` | Shared CAP coefficients |
| `bseOpts` | BSE inversion and bound-state search options |
| `omegaSpec`, `etaSpec` | Frequency grid and common displayed spectral broadening |
| `Aref`, `Pref` | Complex reference amplitude and population `abs(Aref).^2` |
| `P_hard`, `P_cap`, `P_bse` | Cell arrays of population curves on `t` |
| `Jref` | Broadened reference spectrum on `omegaSpec` |
| `J_hard`, `J_cap`, `J_bse` | Cell arrays of broadened spectra on `omegaSpec` |
| `NocScan`, `errBSE_dense` | Dense BSE depth scan and sampled maximum population errors |
| `ocTolerances`, `Noc` | Operational tolerances and first qualifying depths |
| `Results` | MATLAB table exported as termination_results.csv |
| `OperationalDepth` | MATLAB table exported as operational_closure_depths.csv |
| `OperationalScan` | MATLAB table exported as operational_closure_scan.csv |
| `Nstress`, `Phard5`, `Pbse5` | Diagnostic depth and hard-termination/BSE populations on `t` |

For each population or spectrum cell array, entry `{j}` corresponds to `Nscan(j)`. CAP entry `{j}` additionally uses `NdScan(j)`. The grids contain 3001 points each under the defaults. The diagnostic variable names retain the suffix `5`; use the saved `Nstress` to identify the actual depth if parameters were changed.

## Figure files

Each basename below has both `.pdf` and `.png` outputs.

| Basename | Content |
| --- | --- |
| `Figure1_time_domain` | Population curves for the reference and main compared configurations |
| `Figure2_effective_spectrum` | Broadened spectra for the same configurations |
| `Figure3_time_error_vs_N` | Maximum sampled population error versus depth |
| `Figure4_operational_closure_depth` | Dense BSE error scan and operational tolerances |
| `Figure5_spectral_error_vs_N` | Relative spectral L2 error versus depth |
| `Diagnostic_BSE_stress` | Population comparison at the diagnostic depth, default `Nstress = 5` |

The PDF files contain vector graphics and the PNG files are exported at 300 dpi. These rendered figures supplement the numerical data rather than replacing the full-curve archive.
