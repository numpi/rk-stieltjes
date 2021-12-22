% Global options
max_steps = 50;

% Available functions
invsqrt = @(z) 1 ./ sqrt(z);

%sqrtphi = @(z) (1 - exp(-z.^(-.5))) ./ sqrt(z);
%sqrtphi = @(z) (1 - exp(-z.^(-.5))); % Density: sin(1 / sqrt(t))/pi

% Test1: 1 / sqrt(z) 10.000, 50.000, 100.000 A = lapl
% Test2: 1 / sqrt(z), fix 50.000, various spectra (uniform, gap, well-cond)

% Other functions to test (100.000?), A = lapl
sqrtphi = @(z) (1 - exp(-sqrt(z))) ./ z; % Density: sin(sqrt(t))/pi/t
z2 = @(z) z.^(-0.2);
z8 = @(z) z.^(-0.8);
logm = @(z) log(1 + z) ./ z;

phi = @(z) (1 - exp(-z)) ./ z;
lwz = @(z) lambertw(z) ./ z ./ sqrt(z);

% Configurations to test
confs = { ...
    % Inverse square root / Laplacian / various dimensions
    struct('name', '1d-50000-phi', 'n', 50000, 'fun', phi, 'matrix', 'lapl'); ...    
    struct('name', '1d-50000-lwz', 'n', 50000, 'fun', lwz, 'matrix', 'lapl'); ...       
};

for j = 1 : length(confs)
    n = confs{j}.n;
    func  = confs{j}.fun;
    matrix = confs{j}.matrix;
    name = confs{j}.name;
    
    fprintf('Running test with the following configuration:\n');
    fprintf(' > Name: %s\n', name);
    fprintf(' > N: %d\n', n);
    fprintf(' > Matrix: %s\n', matrix);
    fprintf(' > Function: ');
    
    if isequal(func, invsqrt)
        fprintf('Inverse square root\n');
    elseif isequal(func, sqrtphi)
        fprintf('Phi function composed with square root\n');
    elseif isequal(func, logm)
        fprintf('Matrix logarithm\n');
    elseif isequal(func, phi)
        fprintf('Phi function\n');
    elseif isequal(func, lwz)
        fprintf('LabertW / Z^(3/2)\n');
    else
        fprintf('\n');
        % error('Unsupported function');
    end
    
    % Generate a random RHS
    u = randn(n, 1);
    u = u / norm(u);    
        
    % Build the (possibly shifted) Laplacian
    switch matrix
        case 'lapl'
            % Diffusion coefficient
            dc = 0.01;

            % Timestep
            dt = 0.1;

            param = dc * dt * (n + 1).^2;            
            
            A = param * spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n);
            l = param * 4 * sin(pi * (1:n)/(2*(n+1))).^2';            
            %A = spdiags(l, 0, n, n);
            % Compute the exact solution
            % x = func(l) .* (u);
            x = idst(func(l) .* dst(u));
        case 'gap'
            n1 = 20;
            n2 = n - n1;
            l = [ 1e-3 * sin(pi * (1:n1)/(2*(n1+1))).^2' ; 1e1 + 1e3 * sin(pi * (1:n2)/(2*(n2+1))).^2' ];
            % l = [ chebpts(10, [1e-3 1e-3]) ; chebpts(n-10, [10 1000]) ];
            A = spdiags(l, 0, n, n);
            x = func(l) .* u;
        case 'uniform'
            l = linspace(1/n, 1, n).';
            A = spdiags(l, 0, n, n);
            x = func(l) .* u;
        case 'well-conditioned'
            l = 4 * sin(pi * (1:n)/(2*(n+1))).^2' + 1e-3;
            A = spdiags(l, 0, n, n);
            x = func(l) .* u;
    end
    
    % Compute a bound for the spectrum, with slightly relaxed a = a/2
    a = min(l); b = max(l);
    
    % Bound in Corollary 3.14
    rho = sqrt(exp(-pi^2/(log(4*b/a))));
  	gammalk = 2.23 + 2/pi * log(4*(1:max_steps)*sqrt(b/a/pi));

    bound = 8 * gammalk .* rho.^(1:max_steps); 
    
    % Compare different methods
    % 0: Extended
    % 1: Optimal Zolotarev poles
    % 2: EDS distributed poles
    
    tE = tic;
    [xE, errE] = fun_extended_1D(func, A, u, max_steps, eps, x);
    tE = toc(tE);
    fprintf('Extended took %2.1f seconds\n', tE);
    
    errZ = zeros(1, max_steps);
    fprintf('Optimal pole runs:    ');
    for jj = 1 : max_steps
        pZ = laplace_poles(a, b, jj);
        tZ = tic;
        xR = fun_rational_1D(func, A, u, pZ);
        tZ = toc(tZ);
        errZ(jj)  = norm(xR - x);
        fprintf('\b\b\b%3d', jj);
    end
    
    fprintf('\nOptimal Zolotarev poles took %2.1f seconds\n', tZ);
    
    tA = tic;
    [xA, errA] = fun_adaptive_rational_1D(func, A, a, b, u, eps, max_steps, x, 'laplace');    
    tA = toc(tA);
    fprintf('Adaptive Zolotarev poles took %2.1f seconds\n', tA);
    
    %tG = tic;
    %[xG, ~, errG] = invsqrtmv(A, u, eps, max_steps, [], x, func);
    %tG = toc(tG);
    %fprintf('Greedy pole selection by Guettel took %2.1f seconds\n', tG);

    r = [ (1:max_steps)', errE', errZ', errA', bound' ];
    
    dlmwrite(sprintf('%s.dat', name), r, '\t');
    
    figure;
    semilogy(1 : max_steps, r(:,2), 'b-', ...
         1 : max_steps, r(:,3), 'r-', ...
         1 : max_steps, r(:,4), 'c-', ...
         1 : max_steps, r(:,5), 'm-');
    legend('extended', 'laplace', 'adaptive', 'bound');
    
    title(name);    
end
