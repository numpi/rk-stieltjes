function x = fun_rational_1D(f, A, u, poles)
%FUN_EXTENDED Evaluate f(A)u via extended Krylov method.

if size(poles, 1) == 1
	polesA = poles;
else
	polesA = poles.';
end

% Construct the rational Krylov subspace
[VA, KA, HA] = rat_krylov(A, u, [polesA inf]);

Al = HA(1:end-1,:) / KA(1:end-1,:);

% Make them symmetric
% Al=(Al+Al')/2;

% Project the RHS
ul = VA(:,1:end-1)' * u;

% Evaluate the small function
Y = fun_diag_1D(f, Al, ul);

x = VA(:,1:end-1) * Y;

end

