function [g,omega,kappa] = woods_chain_coeffs(Nmax,omega_c,alpha)
%WOODS_CHAIN_COEFFS Analytic Woods particle-chain coefficients for s=1.

n = (1:Nmax).';

g = omega_c*sqrt(alpha);

omega = (omega_c/2) .* ...
    (1 + 1 ./ ((2*n-1).*(2*n+1)));

kappa = omega_c .* sqrt(n.*(n+1)) ./ ...
    (2*(2*n+1));
end
