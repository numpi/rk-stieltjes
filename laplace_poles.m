function poles = laplace_poles(a, b, n)
% LAPLACE_POLES Return the optimal poles for evaluating a Laplace-Stieltjes 
% function on a Hermitian matrix with spectrum in [a, b]. 
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

poles = zeros(n, 1);
for j = 1 : n
	[~, ~, poles(j)] =  ellipj((2 * (n - j) + 1)/(2 * n) * ellipke(1 - a^2 / b^2) , 1 - a^2/b^2);
end
poles = -b * poles;
