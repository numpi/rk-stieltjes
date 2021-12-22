function varargout = fun_rational_2D(f, A, B, u, v, poles)
%FUN_EXTENDED Evaluate f(I \otimes A - B^T \otimes I) vec(u*v') in matrix form using extended Krylov subspaces.

if size(poles, 1) == 1
	polesA = [poles inf];
	polesB = [-poles inf];
else
	polesA = [poles.' inf];
	polesB = [-poles.' inf];
end

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

