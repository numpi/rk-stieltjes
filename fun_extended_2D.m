function varargout = fun_extended_2D(f, A, B, u, v, l)
%FUN_EXTENDED Evaluate f(I \otimes A - B^T \otimes I) vec(u*v') in matrix
% form using extended Krylov subspaces.
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

polesA = zeros(1, l+1); polesA(2:2:l) = inf; polesA(end) = inf;
polesB = polesA;

% Construct the rational Krylov subspaces
[VA, ~, ~] = rat_krylov(A,   u, polesA);
[VB, ~, ~] = rat_krylov(B.', v, polesB);

VA = VA(:, 1:end-1);
VB = VB(:, 1:end-1);

Al = VA' * A * VA;
Bl = VB' * B.' * VB;

% Make them symmetric
Al=(Al+Al')/2;
Bl=(Bl+Bl')/2;

% Project the RHS
ul = VA' * u;
vl = VB' * v;

% Evaluate the small function
Y = fun_diag(f, Al, Bl.', ul, vl);

if nargout >= 3
	varargout{1} = VA;
	varargout{2} = Y;
	varargout{3} = VB;
elseif nargout == 2
	varargout{1} = VA * Y;
	varargout{2} = VB;
else
	error('Unsupported number of output arguments');
end

end

