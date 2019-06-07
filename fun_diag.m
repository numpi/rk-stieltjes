function Y = fun_diag(f, A, B, u, v)
%FUN_DIAG Evaluate f(I \otimes A - B^T \otimes I) vec(u*v') in matrix form.
%
% If f(z) = 1/z, this function solves the equation AX - XB - UV' = 0;
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

if ~ishermitian(A) || ~ishermitian(B)
	error('FUN_DIAG only supports hermitian matrices');
end

[Q1, D1] = eig(A);
[Q2, D2] = eig(B);

C = f( bsxfun(@minus, diag(D1), diag(D2).') );

uu = (Q1' * u);
vv = (Q2' * v);

Y = Q1 * (C .* (uu * vv') * Q2');

end

