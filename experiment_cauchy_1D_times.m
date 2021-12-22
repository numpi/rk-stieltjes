% Global options
max_steps = 50;

% Available functions
func = @(z) 1 ./ sqrt(z);
n = 100000;
        
% Build the (possibly shifted) Laplacian
A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n);
l = 4 * sin(pi * (1:n)/(2*(n+1))).^2';
    
% Generate a random RHS
u = randn(n, 1);
u = u / norm(u);
    
% Compute the exact solution
x = idst(func(l) .* dst(u));
    
% Compute a bound for the spectrum
a = min(l)/2; b = max(l);
    
% Bound in Corollary 3.14
rho = exp(-pi^2/(log(16*b/a)));
bound = 8 * func(a) * rho.^(1:max_steps); 

tols = [ .1 1e-2 1e-3 1e-4 1e-5 1e-6 ];

% [ tols, tguettel, tadapt, textended ]
r = zeros(length(tols), 7);
r(:,1) = tols.';

nrmx = norm(x);

for j = 1 : length(tols)            
    tg = tic;    
    [~,~,err] = invsqrtmv(A, u, tols(j) * nrmx, 1000, [], x, func);
    r(j,2) = toc(tg);
    r(j,2) = timeit(@() invsqrtmv(A, u, tols(j) * nrmx, 1000, [], x, func));
    r(j,3) = length(err);
    fprintf('tol = %d, tguettel = %2.1f seconds, its = %d\n', tols(j), r(j,2), length(err));
    
    tg = tic;
    xi = arrayfun(@(j) cauchy_poles_1D_eds(a, b, j), 1 : 50);
    [~,~,err] = invsqrtmv(A, u, tols(j) * nrmx, 1000, xi, x, func);
    r(j,4) = toc(tg);
    r(j,4) = timeit(@() invsqrtmv(A, u, tols(j) * nrmx, 1000, xi, x, func));
    r(j,5) = length(err);
    fprintf('tol = %d, tadaptive = %2.1f seconds, its = %d\n', tols(j), r(j,4), length(err));
    
    tg = tic;    
    [~,errE] = fun_extended_1D(func, A, u, 2000, tols(j) * nrmx, x);
    % [~,~,errE] = invsqrtmv2(A, u, tols(j), 1000, kron(ones(1, 1000), [eps inf]), x, func);
    r(j,6) = toc(tg);
    r(j,7) = length(errE);
    fprintf('tol = %d, textended = %2.1f seconds, its = %d\n', tols(j), r(j,6), length(errE));    
end

dlmwrite('cauchy-times.dat', r, '\t');