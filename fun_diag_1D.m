function Y = fun_diag_1D(f, A, u)
%FUN_DIAG Evaluate f(A) u by diagonalization. 
%
% References:
% [1] Rational Krylov for Stieltjes matrix functions: convergence and pole 
%     selection, S. Massei and L. Robol, 2019.

if ~ishermitian(A)
	error('FUN_DIAG only supports hermitian matrices');
end

[Q1, D1] = eig(A, 'vector');

Y = Q1 * diag(f(D1)) * (Q1' * u);

end

