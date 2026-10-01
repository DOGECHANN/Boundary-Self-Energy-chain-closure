function SigmaQ = boundary_self_energy_chain( ...
    p,N,g,omega,kappa,gamma,SigmaAfterN)
%BOUNDARY_SELF_ENERGY_CHAIN Inward Rule-I propagation to the emitter.
%
% SigmaAfterN is Sigma_{N+1 -> N}.

if nargin < 6 || isempty(gamma)
    gamma = zeros(N,1);
end

if nargin < 7
    SigmaAfterN = zeros(size(p));
end

gamma = gamma(:);

if numel(gamma) ~= N
    error('gamma must contain N entries.');
end

Sigma = SigmaAfterN;

for n = N-1:-1:1
    Sigma = kappa(n)^2 ./ ...
        (p + 1i*omega(n+1) + gamma(n+1)/2 + Sigma);
end

SigmaQ = g^2 ./ ...
    (p + 1i*omega(1) + gamma(1)/2 + Sigma);
end
