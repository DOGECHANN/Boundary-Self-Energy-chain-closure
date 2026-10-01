%% CAP coefficient scan for the chain-termination benchmark
%
% This script selects ONE common CAP profile for all retained chain lengths.
% It does not re-tune the CAP separately for each N.
%
% The primary score is the geometric mean of the maximum population errors
% over N = 10, 25, 40, 60 on 0 <= t <= 300.  A geometric mean is used so
% that the scan balances errors across several orders of magnitude instead
% of being dominated only by the shortest chain.
%
% The spectral L2 errors are then evaluated for the selected pair as an
% independent frequency-domain check.

clear; clc; close all;

packageDir = fileparts(mfilename('fullpath'));
addpath(packageDir,'-begin');

%% Benchmark parameters
omega_c = 1.0;
alpha   = 0.02;
epsilon = 0.5;

omegaInf = omega_c/2;
kappaInf = omega_c/4;

Nscan  = [10 25 40 60];
NdScan = [5 12 20 30];

Nref = 300;
tMax = 300;
dt = 0.1;
t = 0:dt:tMax;

Ncoeff = max(Nref,max(Nscan)+1);
[g,omega,kappa] = woods_chain_coeffs(Ncoeff,omega_c,alpha);

%% Long-chain reference
fprintf('Building Nref = %d reference...\n',Nref);
[Aref,~] = finite_chain_amplitude_spectral( ...
    Nref,t,epsilon,g,omega,kappa,zeros(Nref,1));
Pref = abs(Aref).^2;

%% Two-stage grid scan
% Coarse scan
coarseA = 0:0.005:0.05;
coarseB = 0.20:0.05:0.90;
[bestA,bestB,bestScore,coarseTable] = scan_grid( ...
    coarseA,coarseB,Nscan,NdScan,t,Pref,epsilon,g,omega,kappa);

fprintf('\nCoarse best: aCAP = %.6g, bCAP = %.6g, score = %.6e\n', ...
    bestA,bestB,bestScore);

% Fine scan around the coarse minimum.  Always include aCAP = 0.
fineAmin = max(0,bestA-0.010);
fineAmax = bestA+0.010;
fineBmin = max(0,bestB-0.10);
fineBmax = bestB+0.10;

fineA = unique([0, fineAmin:0.001:fineAmax]);
fineB = fineBmin:0.01:fineBmax;

[bestA,bestB,bestScore,fineTable] = scan_grid( ...
    fineA,fineB,Nscan,NdScan,t,Pref,epsilon,g,omega,kappa);

% Round the reported coefficients to two decimal places for a simple,
% reproducible benchmark profile, then re-evaluate the errors.
aCAP = round(bestA,2);
bCAP = round(bestB,2);

fprintf('\nFine best before rounding: aCAP = %.6g, bCAP = %.6g\n',bestA,bestB);
fprintf('Adopted common profile:     aCAP = %.2f, bCAP = %.2f\n',aCAP,bCAP);
fprintf('Peak loss gamma_max = %.6g, gamma_max/kappaInf = %.6g\n', ...
    aCAP+bCAP,(aCAP+bCAP)/kappaInf);

%% Time-domain errors for adopted profile
errCAP = zeros(numel(Nscan),1);
for j = 1:numel(Nscan)
    N = Nscan(j);
    Nd = NdScan(j);
    gamma = cap_profile(N,Nd,aCAP,bCAP);
    [Acap,~] = finite_chain_amplitude_spectral( ...
        N,t,epsilon,g,omega,kappa,gamma);
    Pcap = abs(Acap).^2;
    errCAP(j) = max(abs(Pcap-Pref));
end

%% Frequency-domain check
omegaSpec = linspace(0,omega_c,3001);
etaSpec = 2.5*omega_c/Nref;
pSpec = etaSpec - 1i*omegaSpec;

SigmaQ_ref = boundary_self_energy_chain( ...
    pSpec,Nref,g,omega,kappa,zeros(Nref,1),zeros(size(pSpec)));
Jref = real(SigmaQ_ref);
refNorm = trapz(omegaSpec,abs(Jref).^2);

specL2CAP = zeros(numel(Nscan),1);
for j = 1:numel(Nscan)
    N = Nscan(j);
    Nd = NdScan(j);
    gamma = cap_profile(N,Nd,aCAP,bCAP);
    SigmaQ_cap = boundary_self_energy_chain( ...
        pSpec,N,g,omega,kappa,gamma,zeros(size(pSpec)));
    Jcap = real(SigmaQ_cap);
    specL2CAP(j) = sqrt( ...
        trapz(omegaSpec,abs(Jcap-Jref).^2)/refNorm);
end

CAPResult = table(Nscan(:),NdScan(:),errCAP,specL2CAP, ...
    'VariableNames',{'N','CAPOnsetNd','CAPTimeMaxError','CAPSpectralL2'});

disp(CAPResult);
writetable(CAPResult,'cap_selected_profile_results.csv');

% Save both scan stages so the parameter choice is reproducible.
ScanResults = [coarseTable; fineTable];
writetable(ScanResults,'cap_parameter_scan.csv');

fprintf('Primary geometric-mean time-error score = %.6e\n', ...
    exp(mean(log(max(errCAP,realmin)))));

%% Local helper
function [bestA,bestB,bestScore,T] = scan_grid( ...
    aGrid,bGrid,Nscan,NdScan,t,Pref,epsilon,g,omega,kappa)

nRows = numel(aGrid)*numel(bGrid);
rows = zeros(nRows,3+numel(Nscan));
row = 0;
bestScore = inf;
bestA = NaN;
bestB = NaN;

for ia = 1:numel(aGrid)
    for ib = 1:numel(bGrid)
        aCAP = aGrid(ia);
        bCAP = bGrid(ib);
        errors = zeros(1,numel(Nscan));

        for j = 1:numel(Nscan)
            N = Nscan(j);
            Nd = NdScan(j);
            gamma = cap_profile(N,Nd,aCAP,bCAP);
            [Acap,~] = finite_chain_amplitude_spectral( ...
                N,t,epsilon,g,omega,kappa,gamma);
            Pcap = abs(Acap).^2;
            errors(j) = max(abs(Pcap-Pref));
        end

        score = exp(mean(log(max(errors,realmin))));

        row = row+1;
        rows(row,:) = [aCAP,bCAP,score,errors];

        if score < bestScore
            bestScore = score;
            bestA = aCAP;
            bestB = bCAP;
        end
    end
end

varNames = {'aCAP','bCAP','GeometricMeanTimeError'};
for j = 1:numel(Nscan)
    varNames{end+1} = sprintf('TimeErrorN%d',Nscan(j)); %#ok<AGROW>
end
T = array2table(rows,'VariableNames',varNames);
end
