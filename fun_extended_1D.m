function [x, err] = fun_extended_1D(f, A, u, S, tol, xexact)
%FUN_EXTENDED Evaluate f(A)u via extended Krylov method.

SS = ek_struct(A, true);

% Construct the rational Krylov subspace
[VA, KA, HA, params] = ek_krylov(SS, u);
ul = VA(:,1:end-1)' * u;

err = nan(1,2);

for j = 2 : ceil(S/2)
    Al = HA(1:end-1,:) / KA(1:end-1,:);
    % Al = .5 * (Al + Al');
    
    % Project the RHS
    % ul = VA(:,1:end-1)' * u;
    ul(size(VA, 2)-1) = 0;

    % Evaluate the small function
    Y = fun_diag_1D(f, Al, ul);
    x = VA(:,1:end-1) * Y;
    
    if exist('xexact', 'var')
        if j > 1
            err(2*j-1) = err(2*j-2);
        end
        err(2*j) = norm(x - xexact);
        
        if err(2*j) < tol
            return;
        end
    end
    
    %if mod(j, 2) == 1
    %    newpole = 0;
    %else
    %    newpole = inf;
    %end    
    
    [VA, KA, HA, params] = ek_krylov(VA, KA, HA, params);
end

end

