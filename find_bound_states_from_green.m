function [Eb,Zb,info] = find_bound_states_from_green( ...
    greenFun,omegaL,omegaR,Emax,useAnalyticDerivative)
%FIND_BOUND_STATES_FROM_GREEN Real poles outside the continuum band.
% Set useAnalyticDerivative=true for green_bse: it supplies [G,D,Dp].
% A sign change alone can also be a POLE of D (a zero of G), so every
% candidate must pass a denominator residual and positive-residue check.
% The spectral sum-rule check in bse_amplitude_kquad also detects missed
% weight: a finite root-search grid is not a guarantee for arbitrary baths.

if nargin < 4 || isempty(Emax)
    Emax = max(10,10*(omegaR-omegaL));
end
if nargin < 5
    useAnalyticDerivative = false;
end
band = omegaR-omegaL;
scale = max([abs(omegaL),abs(omegaR),band,realmin]);
if band <= 0 || Emax <= 0
    error('Need omegaR>omegaL and Emax>0.');
end
edgeMin = 1e-12*scale;
nLog = 1600;
brackets = [];
for edge = [omegaL omegaR]
    distance = logspace(log10(edgeMin),log10(Emax+abs(edge)),nLog);
    if edge == omegaL
        Egrid = edge-distance;
    else
        Egrid = edge+distance;
    end
    Dgrid = get_denominator(greenFun,-1i*Egrid,useAnalyticDerivative);
    f = imag(Dgrid);
    for j = 1:numel(Egrid)-1
        if isfinite(f(j)) && isfinite(f(j+1)) && ...
                sign(f(j))*sign(f(j+1)) <= 0
            brackets(end+1,:) = sort(Egrid(j:j+1)); %#ok<AGROW>
        end
    end
end
Eb = [];
Zb = [];
rejected = 0;
options = optimset('TolX',1e-13*scale,'Display','off');
for j = 1:size(brackets,1)
    try
        [root,~,exitflag] = fzero(@(E) imag(get_denominator( ...
            greenFun,-1i*E,useAnalyticDerivative)),brackets(j,:),options);
        D = get_denominator(greenFun,-1i*root,useAnalyticDerivative);
        if exitflag <= 0 || ~(root < omegaL || root > omegaR) || ...
                ~isfinite(D) || abs(D) > 1e-8*scale
            rejected = rejected+1;
            continue;
        end
        if ~isempty(Eb) && any(abs(Eb-root) < 1e-10*scale)
            continue;
        end
        if useAnalyticDerivative
            [~,~,Dp] = greenFun(-1i*root);
        else
            % Compatibility fallback: vary E on the SAME physical sheet.
            gap = min(abs(root-[omegaL omegaR]));
            h = min(1e-5*scale,0.01*gap);
            p0 = -1i*root;
            Dp = (get_denominator(greenFun,p0-2i*h,false) ...
                -8*get_denominator(greenFun,p0-1i*h,false) ...
                +8*get_denominator(greenFun,p0+1i*h,false) ...
                -get_denominator(greenFun,p0+2i*h,false))/(12i*h);
        end
        Z = 1/Dp;
        if ~isfinite(Z) || abs(imag(Z)) > 1e-7*max(1,abs(real(Z))) || ...
                real(Z) <= 0 || real(Z) > 1+1e-7
            rejected = rejected+1;
            continue;
        end
        Eb(end+1) = root; %#ok<AGROW>
        Zb(end+1) = real(Z); %#ok<AGROW>
    catch problem
        rejected = rejected+1;
        % Do not hide unexpected programming errors.
        if ~strcmp(problem.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign') && ...
                ~strcmp(problem.identifier,'MATLAB:fzero:InvalidFunctionSupplied')
            rethrow(problem);
        end
    end
end
[Eb,order] = sort(Eb);
Zb = Zb(order);
info.candidateBrackets = size(brackets,1);
info.rejectedCandidates = rejected;
end

function D = get_denominator(greenFun,p,useAnalyticDerivative)
if useAnalyticDerivative
    [~,D] = greenFun(p);
else
    D = 1./greenFun(p);
end
end
