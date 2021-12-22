% 
clear all
close all
n = 1000;

%alpha = .7;
% func = @(z) 1 ./ (-z);
% func = @(z) z.^(-alpha);
%func = @(z) (exp(z.^(.5))-1) ./ sqrt(z);
%func = @(z) exp(-z);
func = @(z) (1 - exp(-sqrt(z))) ./ sqrt(z);
func = @(z) 1./(z.*(1+z.^2));
func = @(z) 1./(1 + z + z.^2 + z.^3 + z.^4 + z.^5 + z.^6);
func = @(z) (1 - exp(-z)) ./ z;
%func = @(z) 1 ./ (z .* (1 + z.^2));
%func = @(z) log(1 + 1./z);
%func = @(z) 1 ./ sinh(z);


% Diffusion coefficient
dc = 0.01;

% Timestep
dt = 0.1;

param = dc * dt * (n + 1).^2;

A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n) * param; % A = A * A;

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
r = zeros(max_steps, 6);

rho = exp(-pi^2/(2*log(4*b/a)));
bound = 32 * func(1e-6) *rho.^[1:max_steps]; % Bound in Corollary 4.2
for j = 1 : max_steps
	fprintf('Step %d\n', j);
	gammalk = 2.23 + 2/pi * log(4*j*sqrt(b/a/pi));
    [VAE, YE, VBE] = fun_extended_2D(func, A, -A, u, v, j);
    [VAP, YP, VBP] = fun_polynomial_2D(func, A, -A, u, v, j);
    poles = laplace_poles(a, b, j);
    [VAR, YR, VBR] = fun_rational_2D(func, A, -A, u, v, poles);
    poles2 = cauchy_poles_2D(a, b, j);
    [VAR2, YR2, VBR2] = fun_rational_2D(func, A, -A, u, v, poles2);
	bound(j) = bound(j) * gammalk;
    poles3 = arrayfun(@(j) laplace_poles_eds(a, b, j), 1 : j);
    [VAR3, YR3, VBR3] = fun_rational_2D(func, A, -A, u, v, poles3);
    r(j, 1:6) = [ bound(j), norm(X - VAE * YE * VBE'), norm(X - VAP * YP * VBP'),  norm(X - VAR * YR * VBR'),  norm(X - VAR2 * YR2 * VBR2'), norm(X - VAR3 * YR3 * VBR3')   ];
end

s = svd(X, 'econ');
 dlmwrite('data/laplace_stieltjes_2D.dat', [ (1:max_steps).', r, s(2:max_steps + 1)  ], '\t');

semilogy(1 : max_steps, r(:,1), 'b-', ...
         1 : max_steps, r(:,2), 'r-', ...
         1 : max_steps, [ s(2 : size(u,2) : size(u,2) * (max_steps + 1)) ] , 'k--', ...
         1 : max_steps, r(:,3), 'c-', ...
	 1 : max_steps, r(:,4), 'g-', ...
	 1 : max_steps, r(:,5), 'm-', ...
     1 : max_steps, r(:,6), 'c--');
legend('bound', 'extended', 'singvals', 'polynomial', 'laplace', 'cauchy', 'eds');
