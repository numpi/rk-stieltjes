function [x, err] = fun_adaptive_rational_1D(f, A, a, b, u, tol, max_steps, true_sol, kind)
%FUN_EXTENDED Evaluate f(A)u via extended Krylov method.

% Vector used to store the errors, assuming that the true solution has been
% passed in the true_sol parameter
err = [];

% These are required to determined the Moebius transform that maps the set
% [-inf, 0] U [a, b] into [-1, -sqrt(ha)] U [sqrt(ha), 1]. iTc is the
% inverse transform. 
% D = sqrt(b^2 - a * b);
% Tc = @(z) (D + z - b)./(D - z + b);
% iTc = @(z) ((b + D) * z + b - D) ./ (1 + z);
% ha = (-Tc(0))^2;
% 
% % EDS part: we determine the function aa(t) and aa1(t) that define the CDF
% % and the density of the probability distribution for the distribution of
% % the poles. 
% Kp = ellipke(1 - ha);
% 
% if isinf(Kp)
%     Kp = abs(log(sqrt(ha))) + 1.386694159834995;
% end
% 
% a  = @(y) -lellipf(asin(sqrt(-y+1.0)),1.0./sqrt(1.0-ha),5e-16) / Kp + 1;
% a1 = @(y) 1./(2*Kp*sqrt((y-ha).*y.*(1-y)));
% 
% poles = iTc(-eds_get_pole(a, a1, ha, 1));

if ~exist('kind', 'var')
    kind = 'cauchy';
end

switch kind
    case 'cauchy'
        poles_func = @cauchy_poles_1D_eds;
    case 'laplace'
        poles_func = @laplace_poles_eds;
end

poles = poles_func(a, b, 1);

% Construct the rational Krylov subspace
[VA, KA, HA, params] = rat_krylov(A, u, [poles inf]);
ul = VA(:,1:end-1)' * u;

%VA = VA(:, 1:end-1);
for j = 1 : max_steps 
    Al = HA(1:end-1, :) / KA(1:end-1,:);
    
	% Make them symmetric
	% Al = (Al + Al')/2;
    
	% Project the RHS
    % ul = VA(:,1:end-1)' * u;
    
	% Evaluate the small function
	Y = fun_diag_1D(f, Al, ul);

	err(j) = norm(true_sol - VA(:, 1:end-1) * Y);
    
    if err(j) < tol        
        x = VA(:, 1:size(Y, 1)) * Y;
        return;
    end        

    % new_pole = iTc(-eds_get_pole(a, a1, ha, j+1));
    new_pole = poles_func(a, b, j+1);
	[VA, KA, HA] = rat_krylov(A, VA, KA, HA, new_pole, params);
    ul = [ ul ; 0 ];
    
    % Simple trick: we want to the infinity pole at the end of the space, 
    % hence, we perform a swap
    G = planerot(KA(end-1:end,end));
    KA(end-1:end,:) = G * KA(end-1:end,:);
    HA(end-1:end,:) = G * HA(end-1:end,:);
    VA(:,end-1:end) = VA(:,end-1:end) * G';    
    ul(end-1:end) = G * ul(end-1:end);
    G = planerot(HA(end,end:-1:end-1)');
    HA(:,end-1:end) = HA(:,end-1:end) * G;
    KA(:,end-1:end) = KA(:,end-1:end) * G;
end

x = VA(:, 1:size(Y, 1)) * Y;

end


