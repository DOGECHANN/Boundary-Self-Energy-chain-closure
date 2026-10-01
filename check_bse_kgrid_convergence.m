%% Optional BSE k-grid convergence check

clear; clc; close all;

omega_c = 1.0;
alpha   = 0.1;
epsilon = 0.25;

omegaInf = omega_c/2;
kappaInf = omega_c/4;

N = 25;

tMax = 2500;
dt = 0.5;
t = 0:dt:tMax;

[g,omega,kappa] = woods_chain_coeffs(N+1,omega_c,alpha);

q1.eta = 1e-10;
q1.phaseStepMax = pi/8;
q1.NkMin = 8001;
q1.findBoundStates = true;
q1.boundEmax = 10*omega_c;

q2 = q1;
q2.phaseStepMax = pi/12;
q2.NkMin = 12001;

[A1,info1] = bse_amplitude_kquad( ...
    t,N,epsilon,g,omega,kappa, ...
    omegaInf,kappaInf,omega_c,q1);

[A2,info2] = bse_amplitude_kquad( ...
    t,N,epsilon,g,omega,kappa, ...
    omegaInf,kappaInf,omega_c,q2);

assert(abs(info1.totalWeight-1)<1e-7);
assert(abs(info2.totalWeight-1)<1e-7);
fprintf('max amplitude difference = %.6e\n',max(abs(A1-A2)));
assert(max(abs(A1-A2))<1e-7,'Increase quadrature resolution.');
P1 = abs(A1).^2;
P2 = abs(A2).^2;

fprintf('Default Nk = %d\n',info1.Nk);
fprintf('Strict  Nk = %d\n',info2.Nk);
fprintf('max population difference = %.6e\n', ...
    max(abs(P1-P2)));

figure('Color','w');
semilogy(t,max(abs(P1-P2),1e-16),'LineWidth',1.5);
xlabel('t');
ylabel('|P_e^{default}-P_e^{strict}|');
title('BSE k-grid convergence');
grid on; box on;
