function [f,xi,err] = invsqrtmv(A,v,tol,maxit,xi,exact,func)
%INVSQRTMV   Computes the matrix inverse square root times vector.
%   F = INVSQRTMV(A,V) returns an approximation F of sqrtm(A)\V for an
%   N x N matrix A and vector V of length N, without forming the matrix
%   inv(sqrtm(A)) explicitly. The method computes a rational Arnoldi 
%   approximation with automatically selected poles based on a heuristic. 
%   The result should therefore be validated, e.g., by running INVSQRTMV
%   twice and looking at the residual norm : 
%          F = INVSQRTMV(A,V); G = INVSQRTMV(A,F); norm(G - A\V)
%
%   F = INVSQRTMV(A,V,TOL) specifies the tolerance of the method. If TOL is 
%   [] then the default value of 1e-8 is used.
%
%   F = INVSQRTMV(A,V,TOL,MAXIT) specifies the max. number of iterations. 
%   If MAXIT is [] then the default value of min(N,50) is used.
%
%   F = INVSQRTMV(A,V,TOL,MAXIT,XI) allows to specify the pole sequence
%   used for building the rational Krylov space. If XI is [] then a 
%   heuristic pole selection strategy is used.
%
%   [F,XI] = INVSQRTMV(A,V,...) returns a vector XI with the poles used for
%   building the rational Krylov space.
%
%   [F,XI,ERR] = INVSQRTMV(A,V,A,V,TOL,MAXIT,XI,EXACT) returns a vector ERR
%   with the error of the rational Arnoldi approximation for every order. 
%   If EXACT is not provided or [] then ERR contains the estimated error.
%
%   Please note that this is a RESEARCH CODE. It is distributed in the hope 
%   that it will be useful, but WITHOUT ANY WARRANTY; without even the 
%   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% 
%   (C) Stefan Guettel, June 2011

if nargin < 6 || isempty(exact), exact = []; end;
if nargin < 5 || isempty(xi), 
    adaptive = 1;
    xi = []; xi_cand = -logspace(-8,8,1601);
else
    adaptive = 0;
end;
N = length(v);
if nargin < 4 || isempty(maxit), maxit = min(N,50); end;
if nargin < 3 || isempty(tol), tol = 1e-8; end;

if ~exist('func', 'var')
    fm = @(X) inv(sqrtm(full(X)));
else
    fm = @(X) fun_diag_1D(func, X, eye(size(X)));
end
if N == 1, f = fm(A)*v; err = 0; return; end;
solver = @(A,b) A\b;    % can use your favorite solver here
zeta = ones(maxit,1);   % may be optimized
beta = norm(v);
V = zeros(N,maxit+1); 
V(:,1) = v/beta;
H = zeros(maxit+1,maxit);
I = speye(N);
reorth = 0; fmR0 = []; xi = xi(:).';

for n = 1:maxit,
    % polynomial Arnoldi step for computing R = V'*A*V
    Av = A*V(:,n) - zeta(n)*V(:,n);
    w = Av;
    for reo = 0:reorth,
        for k = 1:n,
            h = V(:,k)'*w;
            w = w - V(:,k)*h;
            H(k,n) = H(k,n) + h;
        end;
    end;
    K = H(1:n,1:n)*diag([1./xi(1:n-1),0]) + eye(n);
    R = (H(1:n,1:n) + diag(zeta(1:n))) / K;
    
    if adaptive,
        % compute Ritz values of order n
        theta = eig(R(1:n,1:n));
        % look for minimum of nodal rational function sn
        sn = evalrat(theta,[xi(1:n-1),inf],xi_cand);        
        [dummy,ind] = min(abs(sn)); 
        xi(n) = xi_cand(ind);
    end;
        
    % overwrite rational Arnoldi decomposition with new pole xi(n)
    H(:,n) = 0;
    w = solver(I - A/xi(n),Av);
    for reo = 0:reorth,
        for k = 1:n,
            h = V(:,k)'*w;
            w = w - V(:,k)*h;
            H(k,n) = H(k,n) + h;
        end;
    end;
    H(n+1,n) = norm(w);
    w = w/H(n+1,n);
    V(:,n+1) = w;
    
    % compute Arnoldi approximation and error estimate
    fmR = fm(R)*(beta*eye(n,1));
    if ~isempty(exact),
        f = V(:,1:n)*fmR;
        err(n) = norm(f - exact);
    else
        err(n) = norm([fmR0;0]-fmR);
        fmR0 = fmR;
    end;
    if err(n) < tol,
        f = V(:,1:n)*fmR;
        break;
    end;
    if n == maxit,
        f = V(:,1:n)*fmR;
        disp(['Max. iteration number of ' num2str(maxit) ' reached.']);
    end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sn = evalrat( theta,xi,x )
%EVALRAT Evaluate rational nodal function
%           (x^n + ...)
%   sn(x) = -----------
%            1 + ...
% with zeros theta and poles xi at points x
% in a stable way.
 
n = length(theta);
sn = 0*x + 1;
for j = 1:n,
    % find xi closest to min(abs(sn))
    [dummy,ind] = min(abs(sn));
    xm = x(ind);
    [dummy,ind] = min(abs(xi-xm));
    sn = sn ./ (1-x/xi(ind)+eps);
    xi(ind) = NaN;
    
     % find theta closest to max(abs(sn))
    [dummy,ind] = max(abs(sn));
    tm = x(ind);
    [dummy,ind] = min(abs(theta-tm));
    sn = sn .* (x-theta(ind));
    theta(ind) = NaN;
end;

