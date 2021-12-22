function poles = cauchy_poles_2D(a, b, n)
% Return the optimal poles for evaluating a Cauchy-Stieltjes function on [a, b] on a 2D Hermitian matrix (Kronecker sum)
%

D = sqrt(b^2 - a^2);
ta = (D + a - b)/(D - a + b);

% Find the inverse Moebius map that maps [-1, -l] and [l, 1] into intervals 
% containing[-inf, -a] and [a, b]
C = @(z) ((b + D) * z + b - D) ./ (1 + z);

zz = zolotarev_poles(n, ta, 1);

% Determine zeros and poles
poles = C(zz); 

end
