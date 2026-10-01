function gamma = cap_profile(N,Nd,aCAP,bCAP)
%CAP_PROFILE Polynomial distributed complex absorbing profile.
% gamma_n = Theta(x_n) [aCAP*x_n + bCAP*x_n^2],
% with x_n = (n-Nd)/(N-Nd).

if Nd >= N
    error('Nd must satisfy Nd < N.');
end

n = (1:N).';
x = (n-Nd)./(N-Nd);

gamma = zeros(N,1);
mask = x > 0;
gamma(mask) = aCAP*x(mask) + bCAP*x(mask).^2;
end
