%
% Code replicating the experiment with Cauchy-Stieltjes function (no
% Kronecker sums) in [1]. 
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.


n = 1000;
alpha = .5;
func = @(z) z.^(-alpha);

% Optional parameter for scaling the matrix A
param = 1;

A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n) * param;

l = ((2 - 2 * cos(pi * (1:n)./(n+1)))).' * param;
V = sin((1:n)' * (1:n) ./ (n+1) * pi) * sqrt(2 / (n + 1));

u = randn(n, 1);
u = u / norm(u);
x = V * diag(func(l)) * V' * u; % benchmark result
a = min(l); b = max(l);
max_steps = 50;
r = zeros(max_steps, 4);

% Bound in Corollary 3.14
rho = exp(-pi^2/(log(16*b/a)));
bound = 8 * func(a) * rho.^[0:max_steps-1];

for j = 1 : max_steps
	fprintf('Step %d\n', j);
	xE = fun_extended_1D(func, A, u, j);
	xP = fun_polynomial_1D(func, A, u, j);
	poles = cauchy_poles_1D(a, b, j);
	xR = fun_rational_1D(func, A, u, poles);
	poles2 = laplace_poles(a, b, j);
	xR2 = fun_rational_1D(func, A, u, poles2);
	r(j, 1:5) = [ bound(j), norm(x - xE), norm(x - xP),  norm(x - xR),  norm(x - xR2)  ];
end

% Optional instruction for saving data for the paper
% dlmwrite('data/cauchy_stieltjes_1D.dat', [ (1:max_steps).', r], '\t');

semilogy(1 : max_steps, r(:,1), 'b-', ...
	1 : max_steps, r(:,2), 'r-', ...
	1 : max_steps, r(:,3), 'c-', ...
	1 : max_steps, r(:,4), 'g-', ...
	1 : max_steps, r(:,5), 'm-');
legend('bound', 'extended', 'poly', 'cauchy', 'laplace');

