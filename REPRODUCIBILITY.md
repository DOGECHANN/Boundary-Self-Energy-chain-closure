# Reproducing the chain-closure benchmark

This document describes the supplied source code. It is not a record of a completed MATLAB run. Record the MATLAB release, operating system, source commit, and verification output when producing the public data release.

## 1. Run the solver checks

Open the package folder in MATLAB and run:

```matlab
diary('solver_verification_log.txt')
disp(version)
disp(computer)
verification = verify_bse_solver();
writetable(verification,'solver_verification.csv');
diary off
```

The verification function checks spectral normalization, bound-state handling, the analytic Green-function denominator derivative, and quadrature stability. It compares the BSE amplitude with a long finite chain constructed from the **same closure model**. This tests the solver; the main benchmark separately measures the closure error relative to the original mapped reference chain.

The checks use `(alpha, epsilon) = (0.1, 0.25), (0.02, 0.5), (0.3, 0.25)` and `N = [5, 10, 25, 60]`. Assertions require spectral normalization error below `1e-8`, amplitude error against the same closure below `2e-7`, and amplitude changes under the stricter quadrature settings below `1e-7`. A deliberately omitted bound-state contribution must trigger the spectral-weight check. The function prints its success message only after all assertions pass.

## 2. Run and archive the default benchmark

```matlab
diary('benchmark_run_log.txt')
disp(version)
disp(computer)
run_chain_closure_benchmark_fixed
save('benchmark_data.mat', ...
    'omega_c','alpha','epsilon','g','omegaInf','kappaInf', ...
    'omega','kappa','t','tMax','dt','Nref','Nscan','NdScan', ...
    'aCAP','bCAP','bseOpts','etaSpec','omegaSpec', ...
    'Aref','Pref','P_hard','P_cap','P_bse', ...
    'Jref','J_hard','J_cap','J_bse', ...
    'NocScan','errBSE_dense','ocTolerances','Noc', ...
    'Results','OperationalDepth','OperationalScan', ...
    'Nstress','Phard5','Pbse5');
diary off
```

Run the `save` command before another script clears the workspace. `benchmark_data.mat` contains the plotted population and spectrum curves, their grids, and the associated settings. The figure and CSV outputs are written automatically by the benchmark. Keep these files with the source revision used to generate them.

### Physical and numerical settings

Frequencies and couplings are expressed in units of the cutoff frequency, with `omega_c = 1`; time is expressed in inverse-cutoff units. The bath spectral-density convention is

$$J_0(\nu)=2\pi\alpha\nu,\qquad 0\leq\nu\leq 1.$$

| Setting | Default |
| --- | --- |
| Coupling parameter `alpha` | `0.02` |
| Qubit frequency `epsilon` | `0.5` |
| Qubit-to-chain coupling `g` | `sqrt(alpha)` |
| Homogeneous tail frequency `omegaInf` | `0.5` |
| Homogeneous tail coupling `kappaInf` | `0.25` |
| Time grid | `0:0.1:300`, 3001 samples |
| Reference chain length `Nref` | `300` under the default automatic rule |
| Comparison depths `Nscan` | `[10, 25, 40, 60]` |
| CAP onset indices `NdScan` | `[5, 12, 20, 30]` |
| CAP coefficients `(aCAP, bCAP)` | `(0, 0.57)` |
| Frequency grid `omegaSpec` | 3001 equally spaced points on `[0, 1]` |
| Displayed spectral broadening `etaSpec` | `2.5/Nref = 1/120` |
| BSE inversion regulator `bseOpts.eta` | `1e-10` |
| Minimum quadrature nodes `bseOpts.NkMin` | `8001` |
| Maximum target phase step `bseOpts.phaseStepMax` | `pi/8` |
| Spectral-weight tolerance | `1e-8` |
| Bound-state search | Enabled, `boundEmax = 10` |
| Operational closure scan | `N = 1:60` |
| Operational error tolerances | `[1e-3, 1e-4, 1e-5]` |

The BSE inversion regulator and the displayed spectral broadening serve different purposes and have different values.

### Reference chain

Both time-domain and spectrum comparisons use the long finite chain. The automatic reference length is

```matlab
Nref = max(300, ceil(0.5*vMax*1.6*tMax) + 20);
```

Here `vMax = 2*kappaInf = 0.5`. The default gives `Nref = 300` and an estimated endpoint round-trip time of `1200`, beyond the observation window `tMax = 300`. The return-time estimate guides the reference choice; it is not an independent convergence proof.

To fix the reference length explicitly, set `autoNref = false` and `NrefManual = 300` in the main script. The supplied `NrefManual = 1000` is inactive while `autoNref = true`. Changing `Nref` also changes `etaSpec` unless the spectral-broadening rule is edited.

### Closure and CAP conventions

For hard termination and CAP, `N` counts bath sites; the qubit adds one degree of freedom. For the matched BSE closure, `N` labels the closure depth: the code uses the actual bath frequencies through site `N+1` and the actual links `kappa(1:N)`, followed by a homogeneous continuation beyond site `N+1`. This is the convention in `green_bse.m` and in the main script's boundary self-energy.

For a CAP onset `Nd`, define `x_n = (n-Nd)/(N-Nd)`. The absorption rate is zero for `n <= Nd` and is `gamma_n = aCAP*x_n + bCAP*x_n^2` for `n > Nd`. The effective chain Hamiltonian contains `-i*gamma_n/2` on its diagonal. The same CAP coefficients are used for all compared depths and for the spectral comparison.

### Numerical methods and errors

Finite-chain amplitudes are evaluated by eigendecomposition. For the BSE closure, the continuum contribution is integrated on the wave-number grid with composite Simpson quadrature, using `omega(k) = omegaInf - 2*kappaInf*cos(k)`. Bound-state residues are added explicitly. The quadrature node count grows with the observation time and may be refined further if the spectral-weight check fails. The solver checks total spectral weight without renormalizing it to one.

The time-domain error is the largest absolute population difference **on the sampled time grid**:

$$E_t=\max_j\left|P(t_j)-P_{\mathrm{ref}}(t_j)\right|.$$

The spectrum is the real part of the qubit-port BSE evaluated at `p = etaSpec - i*omega`. All schemes and the reference use the same broadening. The reported relative spectral error is

$$E_J=\left[\frac{\int_0^1|J^{(\eta)}(\omega)-J_{\mathrm{ref}}^{(\eta)}(\omega)|^2\,d\omega}{\int_0^1|J_{\mathrm{ref}}^{(\eta)}(\omega)|^2\,d\omega}\right]^{1/2},$$

with both integrals evaluated using `trapz` on `omegaSpec`.

For each tolerance, the operational closure depth is the first qualifying integer in `1:60` whose sampled population error is at most that tolerance. A missing qualifying depth is represented by `NaN` in MATLAB.

## 3. Optional CAP parameter selection

```matlab
scan_cap_coefficients
```

This script uses the same physical parameters, a 300-site reference, and the same time grid as the default benchmark. It minimizes the geometric mean of the four maximum population errors at `N = [10, 25, 40, 60]`.

The coarse grid is `aCAP = 0:0.005:0.05` and `bCAP = 0.20:0.05:0.90`. A fine grid is then centered on the best coarse pair, with steps `0.001` and `0.01`, respectively. The selected coefficients are rounded to two decimal places and reevaluated. The spectrum is evaluated afterward using the selected profile; it is not the optimization objective.

The script exports `cap_parameter_scan.csv` and `cap_selected_profile_results.csv`. It does not change the hard-coded coefficients in the main benchmark.

## 4. Optional long-time quadrature check

```matlab
check_bse_kgrid_convergence
```

This separate check uses `alpha = 0.1`, `epsilon = 0.25`, `N = 25`, and `t = 0:0.5:2500`. It compares two quadrature settings with phase-step targets `pi/8` and `pi/12`, and minimum node counts `8001` and `12001`. The actual node counts also depend on the longer time window. Assertions require both spectral-weight errors and the maximum amplitude difference to be below `1e-7`. The script displays a plot but does not automatically export it.

## Release record

Before publishing results, retain the MATLAB release and platform, source commit or release tag, run logs, generated summary CSVs, and `benchmark_data.mat`. Record any parameter changes alongside their outputs. The documentation and source inspection alone do not establish that a numerical run passed the verification checks.
