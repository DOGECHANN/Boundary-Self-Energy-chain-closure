%% Chain termination benchmark
%
% Main comparison:
%   1) long finite chain as the time-domain reference;
%   2) hard termination;
%   3) CAP optimized for time-domain population error, with chain-dependent onset;
%   4) BSE closure.
%
% The BSE inverse transform uses a k-space quadrature whose node count grows
% automatically with tMax.
%
% The frequency-domain comparison is also chain-based:
%
%   J_eff^(eta)(omega) = Re Sigma_{1->0}(eta-i omega).
%
% No analytic star-continuum solution is used as the main reference.
% Figure 5 reports the relative L2 error of this broadened spectrum against
% the same broadened long-chain reference on 0 <= omega <= omega_c.
% The CAP coefficients are unchanged in the frequency-domain comparison.

clear; clc; close all;
% Run this complete package from its own folder to avoid older helpers.
packageDir = fileparts(mfilename('fullpath'));
addpath(packageDir,'-begin');
clear sigma_uniform_tail green_bse find_bound_states_from_green bse_amplitude_kquad;
fprintf('BSE solver: %s\n',which('bse_amplitude_kquad'));

%% Publication figure formatting
% Match the visual style used in our previous paper: clean standalone panels,
% no in-axis titles, moderate line weights, light grids, and unobtrusive legends.
set(groot,'defaultAxesFontSize',15);
set(groot,'defaultTextFontSize',15);
set(groot,'defaultLegendFontSize',13);
set(groot,'defaultAxesLineWidth',0.8);
set(groot,'defaultAxesLabelFontSizeMultiplier',1.05);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.0);
set(groot,'defaultLineLineWidth',1.5);

bseColor = [0.4660 0.6740 0.1880];
hardColor = [0 0.4470 0.7410];
capColor = [0.8500 0.3250 0.0980];

%% 1. Physical benchmark

% % Fluctuation
% omega_c = 1.0;
% alpha   = 0.1;
% epsilon = 0.25;

% Exponential Decay
omega_c = 1.0;
alpha   = 0.02;
epsilon = 0.5;

omegaInf = omega_c/2;
kappaInf = omega_c/4;
vMax = 2*kappaInf;

%% 2. User-controlled time grid
tMax = 300;          % freely change, e.g. 500 or 2500
dt   = 0.1;
t    = 0:dt:tMax;

%% 3. Retained chain lengths
% Convergence scan used in Fig. 3.
Nscan  = [10 25 40 60];
NdScan = [5 12 20 30];

% Main unequal-depth comparison used in Figs. 1 and 2:
% hard termination and the fine CAP use N = 60, a coarse CAP uses N = 10,
% and the BSE closure uses N = 10.
NhardMain      = 60;
NcapMain       = 60;
NdCapMain      = 30;
NcapCoarseMain = 10;
NdCapCoarseMain = 5;
NbseMain       = 10;

% Retained only as a diagnostic BSE stress test below.
Nstress = 5;

if numel(Nscan) ~= numel(NdScan)
    error('Nscan and NdScan must contain the same number of entries.');
end
if any(NdScan < 1) || any(NdScan >= Nscan)
    error('Each CAP onset must satisfy 1 <= Nd < N.');
end

%% 4. Long-chain reference, independent of Nscan
autoNref = true;

NrefManual = 1000;
NrefMinimum = 300;
referenceSafetyFactor = 1.6;
referenceExtraSites = 20;

if autoNref
    Nref = max(NrefMinimum, ...
        ceil(0.5*vMax*referenceSafetyFactor*tMax) + ...
        referenceExtraSites);
else
    Nref = NrefManual;
end

tReturnRef = 2*Nref/vMax;

if tReturnRef <= tMax
    warning(['Nref endpoint return estimate lies inside the requested ' ...
             'time window. Increase Nref.']);
end

%% 5. CAP profile
% One common CAP profile is used for every retained chain length.
% The coefficients below were selected by scan_cap_coefficients.m for the
% present bath and observation window; they are not re-tuned for each N
% or for the frequency-domain comparison in Figure 5.
% The scan minimizes the geometric mean of the maximum population errors
% over N = 10, 25, 40, 60, and gives an essentially pure quadratic CAP.
aCAP = 0.00*omega_c;
bCAP = 0.57*omega_c;

% The CAP onset is paired with each retained chain length through NdScan.
% (N,Nd) = (10,5), (25,12), (40,20), (60,30).

%% 6. BSE time inversion controls
bseOpts.eta = 1e-10*omega_c;
bseOpts.normalizationTolerance = 1e-8;
bseOpts.phaseStepMax = pi/8;
bseOpts.NkMin = 8001;
bseOpts.findBoundStates = true;
bseOpts.boundEmax = 10*omega_c;

%% 7. Frequency-domain comparison controls
omegaSpec = linspace(0,omega_c,3001);

% A common broadening smooths the discrete long-chain reference spectrum.
specBroadeningFactor = 2.5;
etaSpec = specBroadeningFactor*omega_c/Nref;

%% 8. Generate all required mapped-chain coefficients
NocScan = 1:60;
Ncoeff = max([Nref,max(Nscan)+1,Nstress+1,max(NocScan)+1]);
[g,omega,kappa] = woods_chain_coeffs( ...
    Ncoeff,omega_c,alpha);

fprintf('Benchmark parameters:\n');
fprintf('  omega_c = %.6g, alpha = %.6g, epsilon = %.6g\n', ...
    omega_c,alpha,epsilon);
fprintf('  g = %.8f, omegaInf = %.8f, kappaInf = %.8f\n', ...
    g,omegaInf,kappaInf);
fprintf('  tMax = %.6g, dt = %.6g\n',tMax,dt);
fprintf('  Nref = %d, estimated return time = %.6g\n', ...
    Nref,tReturnRef);
fprintf('  CAP: aCAP = %.6g, bCAP = %.6g, gammaMax = %.6g\n', ...
    aCAP,bCAP,aCAP+bCAP);
fprintf('  CAP gammaMax/kappaInf = %.6g\n', (aCAP+bCAP)/kappaInf);
fprintf('  CAP pairs (N,Nd): ');
for j = 1:numel(Nscan)
    fprintf('(%d,%d) ',Nscan(j),NdScan(j));
end
fprintf('\n');
fprintf(['  Main comparison: hard N = %d, coarse CAP (N,Nd) = (%d,%d), ' ...
         'fine CAP (N,Nd) = (%d,%d), BSE N = %d\n'], ...
    NhardMain,NcapCoarseMain,NdCapCoarseMain,NcapMain,NdCapMain,NbseMain);
fprintf('  spectral broadening eta = %.6g\n\n',etaSpec);

%% 9. Long-chain time-domain reference
fprintf('Building long-chain reference...\n');

[Aref,~] = finite_chain_amplitude_spectral( ...
    Nref,t,epsilon,g,omega,kappa,zeros(Nref,1));

Pref = abs(Aref).^2;

fprintf('  reference diagonalization complete\n\n');

%% 10. Time-domain termination comparison
nCases = numel(Nscan);

P_hard = cell(nCases,1);
P_cap  = cell(nCases,1);
P_bse  = cell(nCases,1);

errHard = zeros(nCases,1);
errCAP  = zeros(nCases,1);
errBSE  = zeros(nCases,1);

bseNk = zeros(nCases,1);
bseWeight = zeros(nCases,1);
bseInitialError = zeros(nCases,1);
bseBoundWeight = zeros(nCases,1);

for j = 1:nCases
    N  = Nscan(j);
    Nd = NdScan(j);

    fprintf('N = %d, CAP onset Nd = %d\n',N,Nd);

    Ahard = finite_chain_amplitude_spectral( ...
        N,t,epsilon,g,omega,kappa,zeros(N,1));

    P_hard{j} = abs(Ahard).^2;
    errHard(j) = max(abs(P_hard{j}-Pref));

    gamma = cap_profile(N,Nd,aCAP,bCAP);

    Acap = finite_chain_amplitude_spectral( ...
        N,t,epsilon,g,omega,kappa,gamma);

    P_cap{j} = abs(Acap).^2;
    errCAP(j) = max(abs(P_cap{j}-Pref));

    [Abse,bInfo] = bse_amplitude_kquad( ...
        t,N,epsilon,g,omega,kappa, ...
        omegaInf,kappaInf,omega_c,bseOpts);

    P_bse{j} = abs(Abse).^2;
    errBSE(j) = max(abs(P_bse{j}-Pref));

    bseNk(j) = bInfo.Nk;
    bseWeight(j) = bInfo.totalWeight;
    bseInitialError(j) = bInfo.normalizationError;
    bseBoundWeight(j) = sum(bInfo.boundWeights);
    fprintf('  Bound energies: '); fprintf('%.12g ',bInfo.boundEnergies); fprintf('\n');
    fprintf('  Bound weights:  '); fprintf('%.12g ',bInfo.boundWeights); fprintf('\n');

    fprintf('  hard max error = %.6e\n',errHard(j));
    fprintf('  CAP  max error = %.6e\n',errCAP(j));
    fprintf(['  BSE  max error = %.6e, Nk = %d, ' ...
             'spectral weight = %.12f\n\n'], ...
        errBSE(j),bInfo.Nk,bInfo.totalWeight);
end

%% 11. Chain-based effective spectral densities
pSpec = etaSpec - 1i*omegaSpec;

SigmaQ_ref = boundary_self_energy_chain( ...
    pSpec,Nref,g,omega,kappa, ...
    zeros(Nref,1),zeros(size(pSpec)));

Jref = real(SigmaQ_ref);

J_hard = cell(nCases,1);
J_cap  = cell(nCases,1);
J_bse  = cell(nCases,1);

specL2Hard = zeros(nCases,1);
specL2CAP  = zeros(nCases,1);
specL2BSE  = zeros(nCases,1);

% Relative spectral error epsilon_X^(J)(N;eta):
% sqrt(integral_0^omega_c |J_X^(eta)-J_ref^(eta)|^2 d omega /
%      integral_0^omega_c |J_ref^(eta)|^2 d omega).
% All spectra use the same etaSpec, frequency interval, and reference.
refNorm = trapz(omegaSpec,abs(Jref).^2);

for j = 1:nCases
    N  = Nscan(j);
    Nd = NdScan(j);

    SigmaQ_hard = boundary_self_energy_chain( ...
        pSpec,N,g,omega,kappa, ...
        zeros(N,1),zeros(size(pSpec)));

    J_hard{j} = real(SigmaQ_hard);

    gamma = cap_profile(N,Nd,aCAP,bCAP);

    SigmaQ_cap = boundary_self_energy_chain( ...
        pSpec,N,g,omega,kappa, ...
        gamma,zeros(size(pSpec)));

    J_cap{j} = real(SigmaQ_cap);

    SigmaInf = sigma_uniform_tail( ...
        pSpec,omegaInf,kappaInf);

    SigmaAfterN = kappa(N)^2 ./ ...
        (pSpec + 1i*omega(N+1) + SigmaInf);

    SigmaQ_bse = boundary_self_energy_chain( ...
        pSpec,N,g,omega,kappa, ...
        zeros(N,1),SigmaAfterN);

    J_bse{j} = real(SigmaQ_bse);

    specL2Hard(j) = sqrt( ...
        trapz(omegaSpec,abs(J_hard{j}-Jref).^2) / refNorm);

    specL2CAP(j) = sqrt( ...
        trapz(omegaSpec,abs(J_cap{j}-Jref).^2) / refNorm);

    specL2BSE(j) = sqrt( ...
        trapz(omegaSpec,abs(J_bse{j}-Jref).^2) / refNorm);
end


%% Figure/caption style note
% Each figure is exported as a separate standalone panel (not combined into
% a MATLAB subplot window).  The manuscript caption should carry the
% descriptive title, as in our previous paper, so no title is drawn inside
% the axes.

%% 12. Figure 1: population comparison
% Deliberately unequal retained depths:
% hard N=60, coarse CAP (N,Nd)=(10,5), fine CAP (N,Nd)=(60,30), BSE N=10.
idxHardMain      = find(Nscan==NhardMain,1);
idxCapMain       = find(Nscan==NcapMain,1);
idxCapCoarseMain = find(Nscan==NcapCoarseMain,1);
idxBSEMain       = find(Nscan==NbseMain,1);

if isempty(idxHardMain) || isempty(idxCapMain) || isempty(idxCapCoarseMain) || isempty(idxBSEMain)
    error('Main-comparison chain lengths must be included in Nscan.');
end
if NdScan(idxCapMain) ~= NdCapMain
    error('NdCapMain does not match NdScan for NcapMain.');
end
if NdScan(idxCapCoarseMain) ~= NdCapCoarseMain
    error('NdCapCoarseMain does not match NdScan for NcapCoarseMain.');
end

fig1 = figure('Color','w','Name','Time-domain termination comparison');
set(fig1,'Position',[100 100 660 470]);

plot(t,Pref,'LineWidth',1.8); hold on;
plot(t,P_hard{idxHardMain},'--','LineWidth',1.4);
plot(t,P_cap{idxCapCoarseMain},'-.','LineWidth',1.4);
plot(t,P_cap{idxCapMain},'-','LineWidth',1.4);
plot(t,P_bse{idxBSEMain},':','Color',bseColor,'LineWidth',1.8);

xlabel('t');
ylabel('P_e(t)');

lgd = legend(sprintf('Long-chain reference, N_{ref}=%d',Nref), ...
       sprintf('Hard cutoff, N=%d',NhardMain), ...
       sprintf('CAP, (N,N_d)=(%d,%d)',NcapCoarseMain,NdCapCoarseMain), ...
       sprintf('CAP, (N,N_d)=(%d,%d)',NcapMain,NdCapMain), ...
       sprintf('BSE closure, N=%d',NbseMain), ...
       'Location','best');
lgd.Box = 'on';
lgd.Color = 'white';

grid on; box on;
ax = gca;
ax.GridAlpha = 0.22;
ax.MinorGridAlpha = 0.10;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xlim([0 tMax]);
ylim([0 1.02]);



exportgraphics(fig1,'Figure1_time_domain.pdf','ContentType','vector');
exportgraphics(fig1,'Figure1_time_domain.png','Resolution',300);

%% 13. Figure 2: effective spectral density
% Use the same unequal-depth comparison as Fig. 1, including the coarse CAP.
fig2 = figure('Color','w','Name','Effective spectral density');
set(fig2,'Position',[100 100 660 470]);

plot(omegaSpec,Jref,'k','LineWidth',1.8); hold on;
plot(omegaSpec,J_hard{idxHardMain},'--','LineWidth',1.4);
plot(omegaSpec,J_cap{idxCapCoarseMain},'-.','LineWidth',1.4);
plot(omegaSpec,J_cap{idxCapMain},'-','LineWidth',1.4);
plot(omegaSpec,J_bse{idxBSEMain},':','Color',bseColor,'LineWidth',1.8);

xlabel('\omega');
ylabel('$J_{\mathrm{eff}}^{(\eta)}(\omega)$','Interpreter','latex');

lgd = legend(sprintf('Long-chain reference, N_{ref}=%d',Nref), ...
       sprintf('Hard cutoff, N=%d',NhardMain), ...
       sprintf('CAP, (N,N_d)=(%d,%d)',NcapCoarseMain,NdCapCoarseMain), ...
       sprintf('CAP, (N,N_d)=(%d,%d)',NcapMain,NdCapMain), ...
       sprintf('BSE closure, N=%d',NbseMain), ...
       'Location','best');
lgd.Box = 'on';
lgd.Color = 'white';

grid on; box on;
ax = gca;
ax.GridAlpha = 0.22;
ax.MinorGridAlpha = 0.10;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xlim([0 omega_c]);

exportgraphics(fig2,'Figure2_effective_spectrum.pdf','ContentType','vector');
exportgraphics(fig2,'Figure2_effective_spectrum.png','Resolution',300);

%% 14. Figure 3: time-domain error versus N
fig3 = figure('Color','w','Name','Time-domain error versus N');
set(fig3,'Position',[100 100 660 470]);

semilogy(Nscan,errHard,'-o','Color',hardColor,'LineWidth',1.5,'MarkerSize',6); hold on;
semilogy(Nscan,errCAP,'-s','Color',capColor,'LineWidth',1.5,'MarkerSize',6);
semilogy(Nscan,errBSE,'-^', 'Color', bseColor, 'LineWidth',1.5,'MarkerSize',6);

xlabel('Chain length $N$','Interpreter','latex');
ylabel('Population error $\varepsilon_X$','Interpreter','latex');

lgd = legend('Hard cutoff', ...
       'CAP (time-domain optimized)', ...
       'BSE closure', ...
       'Location','best');
lgd.Box = 'on';
lgd.Color = 'white';

grid off;
box on;

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

xticks(Nscan);

exportgraphics(fig3,'Figure3_time_error_vs_N.pdf','ContentType','vector');
exportgraphics(fig3,'Figure3_time_error_vs_N.png','Resolution',300);

%% 15. Figure 5: relative spectral error versus N
% Figures 1-4 retain their existing numbering and output filenames.
% If the benchmark arrays are already in the workspace, run only this
% section to redraw/export Figure 5 without repeating the dynamics.
% Suggested caption:
% Relative L2 error of the broadened effective spectral density versus the
% explicit chain length N, using the long-chain reference Nref=300 and the
% common broadening eta=2.5/Nref on 0 <= omega <= 1. The CAP coefficients
% are those optimized for the time-domain population error and are used
% unchanged here. The same CAP onset sites as in Figure 3 are used.
hardColor = [0 0.4470 0.7410];
capColor = [0.8500 0.3250 0.0980];
bseColor = [0.4660 0.6740 0.1880];
figSpec = figure('Color','w','Name','Relative spectral error versus N');
set(figSpec,'Position',[100 100 660 470]);

semilogy(Nscan,specL2Hard,'-o','Color',hardColor,'LineWidth',1.5,'MarkerSize',6,'DisplayName','Hard cutoff'); hold on;
semilogy(Nscan,specL2CAP,'-s','Color',capColor,'LineWidth',1.5,'MarkerSize',6,'DisplayName','CAP (time-domain optimized)');
semilogy(Nscan,specL2BSE,'-^','Color',bseColor,'LineWidth',1.5,'MarkerSize',6,'DisplayName','BSE closure');

xlabel('Chain length $N$','Interpreter','latex');
ylabel('Relative spectral error $\varepsilon_X^{(J)}$','Interpreter','latex');

lgd = legend('Location','best');
lgd.Box = 'on';
lgd.Color = 'white';

grid off; box on;
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xticks(Nscan);

exportgraphics(figSpec,'Figure5_spectral_error_vs_N.pdf','ContentType','vector');
exportgraphics(figSpec,'Figure5_spectral_error_vs_N.png','Resolution',300);

%% 16. Diagnostic: N=5 BSE stress test
Ahard5 = finite_chain_amplitude_spectral( ...
    Nstress,t,epsilon,g,omega,kappa,zeros(Nstress,1));

Phard5 = abs(Ahard5).^2;

[Abse5,~] = bse_amplitude_kquad( ...
    t,Nstress,epsilon,g,omega,kappa, ...
    omegaInf,kappaInf,omega_c,bseOpts);

Pbse5 = abs(Abse5).^2;

figStress = figure('Color','w','Name','BSE stress test');
set(figStress,'Position',[100 100 660 470]);

plot(t,Pref,'LineWidth',1.8); hold on;
plot(t,Phard5,'--','LineWidth',1.4);
plot(t,Pbse5,':','Color',bseColor,'LineWidth',1.8);

xlabel('t');
ylabel('P_e(t)');

legend(sprintf('Long-chain reference, N_{ref}=%d',Nref), ...
       'Hard cutoff', ...
       'BSE closure', ...
       'Location','best','Box','off');

grid on; box on;
ax = gca;
ax.GridAlpha = 0.22;
ax.MinorGridAlpha = 0.10;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xlim([0 tMax]);
ylim([0 1.02]);

xr = xline(4*Nstress,':', ...
    sprintf('t_{refl} \\approx %d',4*Nstress), ...
    'LabelVerticalAlignment','bottom');
xr.HandleVisibility = 'off';
xr.FontSize = 12;

exportgraphics(figStress,'Diagnostic_BSE_stress.pdf','ContentType','vector');
exportgraphics(figStress,'Diagnostic_BSE_stress.png','Resolution',300);

%% 17. Figure 4: operational closure depth
%
% Dense BSE-only scan.  Figure 3 is intentionally left unchanged as the
% three-method comparison.  Here we evaluate every retained depth and mark
% the first depth that reaches each prescribed population-error tolerance.

% Extend the scan to capture the 1e-5 tolerance at alpha=0.1, epsilon=0.25.
% NocScan is defined before coefficient generation in Section 8.
ocTolerances = [1e-3 1e-4 1e-5];

errBSE_dense = zeros(size(NocScan));

fprintf('\nOperational closure-depth scan:\n');

for j = 1:numel(NocScan)
    N = NocScan(j);

    [AbseDense,~] = bse_amplitude_kquad( ...
        t,N,epsilon,g,omega,kappa, ...
        omegaInf,kappaInf,omega_c,bseOpts);

    PbseDense = abs(AbseDense).^2;
    errBSE_dense(j) = max(abs(PbseDense-Pref));

    fprintf('  N = %2d, BSE max error = %.6e\n', ...
        N,errBSE_dense(j));
end

Noc = nan(size(ocTolerances));

for q = 1:numel(ocTolerances)
    idx = find(errBSE_dense <= ocTolerances(q),1,'first');

    if ~isempty(idx)
        Noc(q) = NocScan(idx);
    end
end

fig4 = figure('Color','w','Name','Operational closure depth');
set(fig4,'Position',[100 100 660 470]);

semilogy(NocScan,errBSE_dense,'-o', ...
    'Color',[0 0.4470 0.7410], ...
    'LineWidth',1.5,'MarkerSize',5, ...
    'DisplayName','Boundary-closure population error');

% semilogy(NocScan,errBSE_dense,'-o','Color',bseColor, ...
%     'LineWidth',1.5,'MarkerSize',5, ...
%     'DisplayName','Boundary-closure population error');
% semilogy(NocScan,errBSE_dense,'-o', ...
%     'LineWidth',1.5,'MarkerSize',5);
hold on;

for q = 1:numel(ocTolerances)
    exponent = round(log10(ocTolerances(q)));

    yl = yline(ocTolerances(q),'--', ...
        sprintf('$\\varepsilon = 10^{%d}$',exponent), ...
        'Interpreter','latex', ...
        'LabelHorizontalAlignment','left', ...
        'LabelVerticalAlignment','bottom');
    yl.HandleVisibility = 'off';

    if ~isnan(Noc(q))
        idx = find(NocScan==Noc(q),1);

        plot(Noc(q),errBSE_dense(idx),'s', ...
            'MarkerSize',8,'LineWidth',1.6, ...
            'HandleVisibility','off');

        text(Noc(q)+0.35,ocTolerances(q)*1.25, ...
            sprintf('N_{oc}=%d',Noc(q)), ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom', ...
            'FontSize',13);
    end
end

xlabel('Chain length $N$','Interpreter','latex');
% ylabel('\epsilon_{\rm BSE}');
ylabel('$\varepsilon$','Interpreter','latex');

lgd = legend('Location','best');
lgd.Box = 'on';
lgd.Color = 'white';

exportgraphics(fig4,'Figure4_operational_closure_depth.pdf', ...
    'ContentType','vector');
exportgraphics(fig4,'Figure4_operational_closure_depth.png', ...
    'Resolution',300);

OperationalDepth = table( ...
    ocTolerances(:),Noc(:), ...
    'VariableNames',{'Tolerance','OperationalClosureDepth'});

disp(OperationalDepth);
writetable(OperationalDepth,'operational_closure_depths.csv');

OperationalScan = table( ...
    NocScan(:),errBSE_dense(:), ...
    'VariableNames',{'N','BSETimeMaxError'});

writetable(OperationalScan,'operational_closure_scan.csv');

%% 18. Numerical summary
Results = table( ...
    Nscan(:),NdScan(:), ...
    errHard,errCAP,errBSE, ...
    specL2Hard,specL2CAP,specL2BSE, ...
    bseNk,bseWeight,bseInitialError,bseBoundWeight, ...
    'VariableNames', ...
    {'N','CAPOnsetNd', ...
     'HardTimeMaxError','CAPTimeMaxError','BSETimeMaxError', ...
     'HardSpectralL2','CAPSpectralL2','BSESpectralL2', ...
     'BSENk','BSESpectralWeight','BSENormalizationError','BSEBoundWeight'});

disp(Results);
writetable(Results,'termination_results.csv');

fprintf('\nLong-chain reference:\n');
fprintf('  Nref = %d\n',Nref);
fprintf('  estimated round-trip return = %.6g\n',tReturnRef);
fprintf('  requested tMax = %.6g\n',tMax);
