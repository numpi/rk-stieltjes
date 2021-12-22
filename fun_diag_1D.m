function Y = fun_diag_1D(f, A, u)
%FUN_DIAG Evaluate f(I \otimes A - B^T \otimes I) vec(u*v') in matrix form.
%
% If f(z) = 1/z, this function solves the equation AX - XB - UV' = 0;

%if ~ishermitian(A) 
%    error('FUN_DIAG only supports hermitian matrices');
%end

[Q1, D1] = eig(A, 'vector');

Y = Q1 * diag(f(D1)) * (Q1 \ u);

end

