function poles = cauchy_poles_1D(a, b, n)
% Return the optimal poles for evaluating a Cauchy-Stieltjes function 
% on [a, b] on a 1D Hermitian matrix (no Kronecker sums). 
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019. 

D = sqrt(b^2 - a * b);
ta = (D + a - b)/(D - a + b);

% Find the inverse Moebius map that maps [-1, -ta] and [ta, 1] into intervals
% containing[-inf, 0] and [a, b]
C = @(z) ((b + D) * z + b - D) ./ (1 + z);

zz = zolotarev_poles(n, ta, 1);

% Determine the poles
poles = C(zz);

