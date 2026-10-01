function results = verify_bse_solver()
%VERIFY_BSE_SOLVER Independent checks for the corrected inversion.
% Does not generate figures. Run before the main benchmark.
% The independent comparator is a long Hermitian chain with EXACTLY
% the same actual prefix and homogeneous continuation as the BSE closure.
% It is a validation comparator, not a replacement BSE algorithm.

omega_c = 1;
omegaInf = 0.5;
kappaInf = 0.25;
t = 0:0.2:300;
L = 300;
[g,omega,kappa] = woods_chain_coeffs(L,omega_c,0.1);
Nvalues = [5 10 25 60];
parameters = [0.1 0.25; 0.02 0.5; 0.3 0.25];
options.NkMin = 8001;
options.eta = 1e-10;
options.findBoundStates = true;
options.normalizationTolerance = 1e-8;
rows = [];

% Branch continuation and analytic derivative above AND below the band.
for E = [-0.1 1.01]
    p0 = -1i*E;
    h = 1e-7;
    [~,~,Dp] = green_bse(p0,10,0.25,g,omega,kappa,omegaInf,kappaInf);
    [~,Dplus] = green_bse(p0+h,10,0.25,g,omega,kappa,omegaInf,kappaInf);
    [~,Dminus] = green_bse(p0-h,10,0.25,g,omega,kappa,omegaInf,kappaInf);
    assert(abs((Dplus-Dminus)/(2*h)-Dp) < 1e-6*max(1,abs(Dp)), ...
        'Physical-sheet continuation or analytic derivative failed.');
end

for q = 1:size(parameters,1)
    alpha = parameters(q,1);
    epsilon = parameters(q,2);
    [g,omega,kappa] = woods_chain_coeffs(L,omega_c,alpha);
    for N = Nvalues
        [Abse,info] = bse_amplitude_kquad(t,N,epsilon,g,omega,kappa, ...
            omegaInf,kappaInf,omega_c,options);
        hybridOmega = omega;
        hybridKappa = kappa;
        hybridOmega(N+2:end) = omegaInf;
        hybridKappa(N+1:end) = kappaInf;
        Ahybrid = finite_chain_amplitude_spectral( ...
            L,t,epsilon,g,hybridOmega,hybridKappa,zeros(L,1));
        amplitudeError = max(abs(Abse-Ahybrid));
        populationError = max(abs(abs(Abse).^2-abs(Ahybrid).^2));
        assert(info.normalizationError < 1e-8,'Spectral normalization failed.');
        assert(amplitudeError < 2e-7,'BSE disagrees with its independent closure.');
        if q == 1
            assert(numel(info.boundEnergies)==1 && info.boundEnergies(1)>1);
        elseif q == 3
            assert(numel(info.boundEnergies)==2, ...
                'Strong-coupling case should include both exterior bound states.');
        end
        rows(end+1,:) = [alpha epsilon N info.totalWeight ...
            numel(info.boundEnergies) amplitudeError populationError]; %#ok<AGROW>
        fprintf('alpha=%.2g epsilon=%.2g N=%d: max |dA|=%.3e, weight=%.12f\n', ...
            alpha,epsilon,N,amplitudeError,info.totalWeight);
    end
end

% A denser grid and smaller broadening separate inversion error from
% finite-depth closure error. This is not a comparison with the exact bath.
[g,omega,kappa] = woods_chain_coeffs(L,omega_c,0.1);
[A1,~] = bse_amplitude_kquad(t,10,0.25,g,omega,kappa, ...
    omegaInf,kappaInf,omega_c,options);
strict = options;
strict.NkMin = 16001;
strict.eta = 1e-11;
[A2,~] = bse_amplitude_kquad(t,10,0.25,g,omega,kappa, ...
    omegaInf,kappaInf,omega_c,strict);
assert(max(abs(A1-A2)) < 1e-7,'BSE quadrature/broadening is not converged.');

% Missing bound-state weight must cause an error, never renormalization.
missing = options;
missing.findBoundStates = false;
expectedFailure = false;
try
    bse_amplitude_kquad(0,10,0.25,g,omega,kappa, ...
        omegaInf,kappaInf,omega_c,missing);
catch problem
    if strcmp(problem.identifier,'BSE:SpectralWeight')
        expectedFailure = true;
    else
        rethrow(problem);
    end
end
assert(expectedFailure,'Missing spectral weight was not detected.');
results = array2table(rows,'VariableNames', ...
    {'alpha','epsilon','N','TotalWeight','BoundStateCount', ...
     'AmplitudeErrorVsSameClosure','PopulationErrorVsSameClosure'});
disp(results);
fprintf('All BSE solver checks passed.\n');
end
