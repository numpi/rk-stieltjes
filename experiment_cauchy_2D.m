% 
clear all
close all
n = 1000;

alpha = .5;
func = @(z) z.^(-alpha);
%func = @(z) (1 - exp(-z.^(-.5))) .* sqrt(z);
%func = @(z) exp(-z);
%func = @(z) (1 - exp(-z)) ./ z;
%func = @(z) 1 ./ (z .* (1 + z.^2));
%func = @(z) log(1 + 1./z);
%func = @(z) 1 ./ sinh(z);
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
r = zeros(max_steps, 7);

rho = exp(-pi^2/(log(8*b/a)));
bound = 4 * (1 + b/a) * func(2 * a) * rho.^[1:max_steps]; % Bound in Corollary 4.9
bound2 = 4 * func(2 * a) * rho.^[1:max_steps]; % Bound in Corollary 4.9

for j = 1 : max_steps
	fprintf('Step %d\n', j);
    [VAE, YE, VBE] = fun_extended_2D(func, A, -A, u, v, j);
    [VAP, YP, VBP] = fun_polynomial_2D(func, A, -A, u, v, j);
    poles = cauchy_poles_2D(a, b, j);
    [VAR, YR, VBR] = fun_rational_2D(func, A, -A, u, v, poles);
    poles2 = laplace_poles(a, b, j);
    [VAR2, YR2, VBR2] = fun_rational_2D(func, A, -A, u, v, poles2);
    poles3 = arrayfun(@(j) cauchy_poles_2D_eds(a, b, j), 1 : j);
    [VAR3, YR3, VBR3] = fun_rational_2D(func, A, -A, u, v, poles3);
	
    r(j, 1:7) = [ bound(j), norm(X - VAE * YE * VBE'), norm(X - VAP * YP * VBP'),  norm(X - VAR * YR * VBR'),  norm(X - VAR2 * YR2 * VBR2'), norm(X - VAR3 * YR3 * VBR3'), bound2(j)];
end

s = svd(X, 'econ');
 dlmwrite('data/cauchy_stieltjes_2D.dat', [ (1:max_steps).', r, s(2:max_steps + 1)  ], '\t');

semilogy(1 : max_steps, r(:,1), 'b-', ...
         1 : max_steps, r(:,2), 'r-', ...
         1 : max_steps, [ s(2 : size(u,2) : size(u,2) * (max_steps + 1)) ] , 'k--', ...
         1 : max_steps, r(:,3), 'c-', ...
         1 : max_steps, r(:,4), 'g-',...
         1 : max_steps, r(:,5), 'm-', ...
         1 : max_steps, r(:,6), 'c--', ...
         1 : max_steps, r(:,7), 'm--');
legend('bound', 'extended', 'singvals', 'polynomial', 'cauchy', 'laplace', 'eds', 'bound2');

