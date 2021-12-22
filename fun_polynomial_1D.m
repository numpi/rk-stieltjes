function x = fun_polynomial_1D(f, A, u, l)
%FUN_EXTENDED Evaluate f(A)u via extended Krylov method.

polesA = inf(1, l+1); 

% Construct the rational Krylov subspace
[VA, ~, ~] = rat_krylov(A,   u, polesA);

VA = VA(:, 1:end-1);
if isstruct(A)
    Al = VA' * A.multiply(1, 0, VA);
else
    Al = VA' * A * VA;
end


% Make them symmetric
Al = (Al+Al')/2;

% Project the RHS
ul = VA' * u;

% Evaluate the small function
Y = fun_diag_1D(f, Al, ul);

x = VA * Y;

end

