function x = fun_polynomial_1D(f, A, u, l)
%FUN_EXTENDED Evaluate f(A)u via a polynomial Krylov method.
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

polesA = inf(1, l+1);

% Construct the rational Krylov subspace
[VA, ~, ~] = rat_krylov(A,   u, polesA);

VA = VA(:, 1:end-1);

Al = VA' * A * VA;


% Make them symmetric
Al = (Al+Al')/2;

% Project the RHS
ul = VA' * u;

% Evaluate the small function
Y = fun_diag_1D(f, Al, ul);

x = VA * Y;

end
