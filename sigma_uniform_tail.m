function [SigmaInf,dSigmaInf] = sigma_uniform_tail(p,omegaInf,kappaInf)
%SIGMA_UNIFORM_TAIL Physical surface response of a uniform half-chain.
% The branch is continued across the imaginary axis OUTSIDE the band.
% On the band, Re(p)=0 denotes the limit from Re(p)>0.
% Optional derivative is analytic away from the two band edges.

if ~isscalar(kappaInf) || ~isreal(kappaInf) || kappaInf <= 0
    error('kappaInf must be a positive real scalar.');
end
z = p + 1i*omegaInf;
rootTerm = sqrt(z.^2 + 4*kappaInf^2);

% For the physical resolvent, rootTerm follows z at large |z|.
% Re(rootTerm)>=0 is correct ONLY in the right half-plane.
left = real(z) < 0;
rootTerm(left) = -rootTerm(left);
onAxis = real(z) == 0;
wrongImag = onAxis & (imag(rootTerm).*imag(z) < 0);
rootTerm(wrongImag) = -rootTerm(wrongImag);

% Rationalization avoids cancellation in (-z+rootTerm)/2.
SigmaInf = 2*kappaInf^2 ./ (z + rootTerm);
if nargout > 1
    dSigmaInf = -SigmaInf ./ rootTerm;
end
end
