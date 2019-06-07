function x = fun_rational_1D(f, A, u, poles)
%FUN_EXTENDED Evaluate f(A)u via rational Krylov method.
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

if size(poles, 1) == 1
	polesA = [poles inf];
else
	polesA = [poles.' inf];
end

% Construct the rational Krylov subspace
[VA, ~, ~] = rat_krylov(A,   u, polesA);

VA = VA(:, 1:end-1);

Al = VA' * A * VA;


% Make them symmetric
Al=(Al+Al')/2;

% Project the RHS
ul = VA' * u;

% Evaluate the small function
Y = fun_diag_1D(f, Al, ul);

x = VA * Y;

end

