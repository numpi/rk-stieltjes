% 
clear all
close all
n = 500;
func = @(z) cos(z);

param = 1;
S = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n) * param;  % A = A * A;
A = inv(S);
l = ((2 - 2 * cos(pi * (1:n)./(n+1)))).' * param; l=l.^(-1);
%l = linspace(.1,2e2, n);
V = sin((1:n)' * (1:n) ./ (n+1) * pi) * sqrt(2 / (n + 1));
%A = V*diag(l)*V';

u = randn(n, 1);
u = u / norm(u);
x = V * diag(func(l)) * V' * u; % benchmark result
a = min(l); b = max(l);
max_steps = 100;
r = zeros(max_steps, 4);

for j = 1 : max_steps
    fprintf('Step %d\n', j);
    xE = fun_extended_1D(func, A, u, j);
    xP = fun_polynomial_1D(func, A, u, j);
    poles = 1i*cauchy_poles_1D(a, b, j);
    poles = zeros(1, j);
    xR = fun_rational_1D(func, A, u, poles);
    poles2 = 1i*laplace_poles(a, b, j);
    %poles2 = 1i*param*rand(1, j);
    xR2 = fun_rational_1D(func, A, u, poles2);
    r(j, 1:4) = [norm(x - xE), norm(x - xP),  norm(x - xR),  norm(x - xR2)  ];

end

dlmwrite('data/cosine.dat', [ (1:max_steps).', r], '\t');

semilogy(1 : max_steps, r(:,1), 'b-', ...
         1 : max_steps, r(:,2), 'r-', ...
         1 : max_steps, r(:,3), 'c-', ...
		 1 : max_steps, r(:,4), 'g-');
legend('extended', 'poly', 'cauchy', 'laplace');

