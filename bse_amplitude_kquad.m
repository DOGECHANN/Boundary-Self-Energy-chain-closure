function [A0,info] = bse_amplitude_kquad( ...
    t,N,epsilon,g,omega,kappa,omegaInf,kappaInf,omega_c,opts)
%BSE_AMPLITUDE_KQUAD Time-domain BSE amplitude using a k-space quadrature.
%
% omega(k)=omegaInf-2*kappaInf*cos(k), 0<=k<=pi.
%
% The node count grows automatically with tMax so that
%
%   tMax*vMax*Delta k <= phaseStepMax.
%
% Phase resolution controls the oscillatory transform, but the continuum
% AND every bound-state residue are needed for a normalized amplitude.
% The sum rule is checked without rescaling the result.

if nargin < 10
    opts = struct();
end

eta = get_opt(opts,'eta',1e-10*omega_c);
normalizationTolerance = get_opt(opts,'normalizationTolerance',1e-8);
maxRefinements = get_opt(opts,'maxRefinements',3);
phaseStepMax = get_opt(opts,'phaseStepMax',pi/8);
NkMin = get_opt(opts,'NkMin',8001);
findBoundStates = get_opt(opts,'findBoundStates',true);
boundEmax = get_opt(opts,'boundEmax',10*omega_c);

t = t(:).';
if isempty(t) || any(~isfinite(t)) || any(t < 0)
    error('t must contain finite nonnegative times.');
end
if eta <= 0 || phaseStepMax <= 0 || NkMin < 3 || ...
        normalizationTolerance <= 0
    error('Invalid BSE quadrature controls.');
end
tMax = max(t);

vMax = 2*kappaInf;

NkPhase = ceil(pi*vMax*max(tMax,1)/phaseStepMax) + 1;
Nk = ceil(max(NkMin,NkPhase));

if mod(Nk,2) == 0
    Nk = Nk+1;
end

k = linspace(0,pi,Nk).';
dk = pi/(Nk-1);

omegaK = omegaInf - 2*kappaInf*cos(k);
jacobian = 2*kappaInf*sin(k);

Gfun = @(p) green_bse( ...
    p,N,epsilon,g,omega,kappa,omegaInf,kappaInf);

G = Gfun(eta-1i*omegaK);
rho = real(G)/pi;

smallNegative = rho < 0 & rho > -1e-11;
rho(smallNegative) = 0;

if any(~isfinite(rho)) || any(rho < -1e-8)
    error('BSE:InvalidDensity','Invalid BSE continuum density.');
end

simpson = ones(Nk,1);
simpson(2:2:end-1) = 4;
simpson(3:2:end-2) = 2;

weight = (dk/3) .* simpson .* jacobian .* rho;

% Locate discrete spectral weight before integrating in time.
Eb = [];
Zb = [];
if findBoundStates
    [Eb,Zb] = find_bound_states_from_green( ...
        Gfun,omegaInf-2*kappaInf, ...
        omegaInf+2*kappaInf,boundEmax,true);
end
info.Nk = Nk;
info.dk = dk;
info.phaseStepAtTmax = tMax*vMax*dk;
info.continuumWeight = sum(weight);
info.boundEnergies = Eb;
info.boundWeights = Zb;
info.totalWeight = info.continuumWeight + sum(Zb);
info.normalizationError = abs(info.totalWeight-1);
info.initialPopulation = abs(info.totalWeight)^2;
info.eta = eta;
info.refinementCount = 0;
if info.normalizationError > normalizationTolerance
    if maxRefinements > 0 && findBoundStates
        % Resolve a narrow continuum feature; do not rescale its weight.
        refined = opts;
        refined.NkMin = 2*(Nk-1)+1;
        refined.maxRefinements = maxRefinements-1;
        [A0,info] = bse_amplitude_kquad(t,N,epsilon,g,omega,kappa, ...
            omegaInf,kappaInf,omega_c,refined);
        info.refinementCount = info.refinementCount+1;
        return;
    end
    error('BSE:SpectralWeight', ...
        ['Spectral weight is %.12g (expected 1). Do not normalize this ' ...
         'away. Check bound poles/search range, then eta and NkMin.'], ...
        info.totalWeight);
end

Acont = zeros(size(t));

uniformGrid = false;

if numel(t) == 1
    uniformGrid = abs(t(1)) < 100*eps;
elseif abs(t(1)) < 100*eps
    dt = t(2)-t(1);
    uniformGrid = max(abs(diff(t)-dt)) < ...
        1e-10*max(1,abs(dt));
end

if uniformGrid
    if numel(t) == 1
        Acont(1) = sum(weight);
    else
        dt = t(2)-t(1);

        phase = ones(Nk,1);
        phaseStep = exp(-1i*omegaK*dt);

        for m = 1:numel(t)
            Acont(m) = weight.'*phase;
            phase = phase.*phaseStep;
        end
    end
else
    chunkSize = 200;

    for i1 = 1:chunkSize:numel(t)
        i2 = min(numel(t),i1+chunkSize-1);
        phase = exp(-1i*(omegaK*t(i1:i2)));
        Acont(i1:i2) = weight.'*phase;
    end
end

Abound = zeros(size(t));
for j = 1:numel(Eb)
    Abound = Abound + Zb(j)*exp(-1i*Eb(j)*t);
end
A0 = Acont + Abound;
end

function value = get_opt(opts,name,default)
if isfield(opts,name)
    value = opts.(name);
else
    value = default;
end
end
