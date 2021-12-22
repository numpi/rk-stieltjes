function pole = laplace_poles_1D_eds(a, b, n)
% Return the optimal poles for evaluating a Laplace-Stieltjes function on [a, b] on a 1D Hermitian matrix (no Kronecker sums)

ta = a/b;
ha = ta^2;

% Find the inverse Moebius map that maps [-1, -ta] and [ta, 1] into intervals 
% containing [-inf, 0] and [a, b]
C = @(z) b * z;

% EDS part: we determine the function aa(t) and aa1(t) that define the CDF
% and the density of the probability distribution for the distribution of
% the poles. 

Kp = ellipke(1 - ha);

if isinf(Kp)
    Kp = abs(log(sqrt(ha))) + 1.386694159834995;
end

a  = @(y) -lellipf(asin(sqrt(-y+1.0)),1.0./sqrt(1.0-ha),5e-16) / Kp + 1;
a1 = @(y) 1./(2*Kp*sqrt((y-ha).*y.*(1-y)));

pole = C(-eds_get_pole(a, a1, ha, n));

end
