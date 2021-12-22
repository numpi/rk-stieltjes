% Code replicating the experiment with Laplace-Stieltjes function (no 
% Kronecker structured matrices) in [1]. 
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

n = 1000;

func = @(z) (1 - exp(-z)) ./ z;
func = @(z) exp(-z)

% Diffusion coefficient
dc = 0.01;

% Timestep
dt = 0.1;

% Scaling of the matrix A
param = dc * dt * (n + 1).^2;

A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n) * param;

l = ((2 - 2 * cos(pi * (1:n)./(n+1)))).' * param;
V = sin((1:n)' * (1:n) ./ (n+1) * pi) * sqrt(2 / (n + 1));


u = randn(n, 1);
u = u / norm(u);
x = V * diag(func(l)) * V' * u; % benchmark result
a = min(l); b = max(l);
max_steps = 50;
r = zeros(max_steps, 4);

% Bound in Corollary 3.13

% Shifting parameter
ss = 5;

rho = exp(-pi^2/(2*log(4*(b + ss) / (a + ss))));
bound = 4 * func(a) * rho.^[0:max_steps-1];

for j = 1 : max_steps
	fprintf('Step %d\n', j);
	xE = fun_extended_1D(func, A, u, j);
	xP = fun_polynomial_1D(func, A, u, j);
	
	poles = laplace_poles(a + ss, b + ss, j);
	xR = fun_rational_1D(@(z) func(z), A, u, poles - ss);
	poles2 = cauchy_poles_1D(a + ss, b + ss, j);
	xR2 = fun_rational_1D(func, A, u, poles2 - ss);
	r(j, 1:5) = [ bound(j), norm(x - xE), norm(x - xP),  norm(x - xR),  norm(x - xR2)  ];
end

% Optional instruction for saving data for the paper
% dlmwrite('data/laplace_stieltjes_1D.dat', [ (1:max_steps).', r], '\t');

semilogy(1 : max_steps, r(:,1), 'b-', ...
	1 : max_steps, r(:,2), 'r-', ...
	1 : max_steps, r(:,3), 'c-', ...
	1 : max_steps, r(:,4), 'g-', ...
	1 : max_steps, r(:,5), 'm-');
legend('bound', 'extended', 'poly', 'laplace', 'cauchy');

