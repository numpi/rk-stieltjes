% Code replicating the experiment with Cauchy-Stieltjes function for
% Kronecker structured matrices in [1]. 
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

n = 1000;

alpha = .5;
func = @(z) z.^(-alpha);

% Optional scaling parameter
param = 1;

A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n) * param;

l = ((2 - 2 * cos(pi * (1:n)./(n+1)))).' * param;
V = sin((1:n)' * (1:n) ./ (n+1) * pi) * sqrt(2 / (n + 1));


u = randn(n, 1);
v = randn(n, 1);


u = u / norm(u);
v = v / norm(v);
% Transform u and v
uu = V' * u;
vv = V' * v;

L = l + l.';
XX = func(L) .* (uu * vv');
X = V * XX * V';
a = min(l); b = max(l);

max_steps = 50;
r = zeros(max_steps, 5);

rho = exp(-pi^2/(log(8*b/a)));
bound = 4 * (1 + b/a) * func(2 * a) * rho.^[1:max_steps]; % Bound in Corollary 4.9
for j = 1 : max_steps
	fprintf('Step %d\n', j);
	[VAE, YE, VBE] = fun_extended_2D(func, A, -A, u, v, j);
	[VAP, YP, VBP] = fun_polynomial_2D(func, A, -A, u, v, j);
	poles = cauchy_poles_2D(a, b, j);
	[VAR, YR, VBR] = fun_rational_2D(func, A, -A, u, v, poles);
	poles2 = laplace_poles(a, b, j);
	[VAR2, YR2, VBR2] = fun_rational_2D(func, A, -A, u, v, poles2);
	
	r(j, 1:5) = [ bound(j), norm(X - VAE * YE * VBE'), norm(X - VAP * YP * VBP'),  norm(X - VAR * YR * VBR'),  norm(X - VAR2 * YR2 * VBR2')];
end

s = svd(X, 'econ');

% Optionally save the data for inclusion in the paper
% dlmwrite('data/cauchy_stieltjes_2D.dat', [ (1:max_steps).', r, s(2:max_steps + 1)  ], '\t');

semilogy(1 : max_steps, r(:,1), 'b-', ...
	1 : max_steps, r(:,2), 'r-', ...
	1 : max_steps, s(2 : size(u,2) : size(u,2) * (max_steps + 1)), 'k--', ...
	1 : max_steps, r(:,3), 'c-', ...
	1 : max_steps, r(:,4), 'g-',...
	1 : max_steps, r(:,5), 'm-');
legend('bound', 'extended', 'singvals', 'polynomial', 'cauchy', 'laplace');

