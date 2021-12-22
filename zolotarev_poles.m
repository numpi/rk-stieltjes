function [p, q] = zolotarev_poles(n, a, b, c, d)
%ZOLOTAREV_POLES Determine poles for the optimal Zolotarev approximant. 
%
%
if ~exist('c', 'var')
    c = a;
    d = b;
end
cx = (a - c) /2; 
if a - cx > 0
	flag_sign = 1;
	a = a - cx;
	b = b - cx;
	c = c + cx;
	d = d + cx;
else
	cx = (-b + d)/2; 
	if -b - cx < 0
		error('Unsupported configuration');
	end
	flag_sign = -1;
	tmp = a;
	a = -b - cx;
	b = -tmp - cx;
	tmp = c;
	c = -d + cx;
	d = -tmp + cx;
end

if n < 1
	rho = exp(-pi^2/(log(4*max(b/a, d/c))));
	n = ceil(log(n/4) / log(rho));
end

eta = 2 * (b - a) * (d - c) / ((a + c) * (b + d));
kp  = 1 / (1 + eta + sqrt(eta * (eta + 2)));

if a == c && b == d
    g = 0.0;
else
    g   = 2 * (kp * (b + d) - (a + c)) / ...
        ((a + c) * (b - d) + kp * (b + d) * (c - a));
end

f = (2 + g * (b - d)) / (b + d);

h = kp * (c - a + 2 * a*c*g) / (a + c);

p = zeros(1, n);
q = zeros(1, n);
k = (1 + kp)*(1 - kp); % 1 - kp^2;

K = ellipk(k);

if isinf(K)
    K = abs(log(kp)) + 1.386694159834995;
end

[~,~,w] = ellipj((2*(1:n)-1)*K/(2*n), k);
p = -(h + w) ./ (w.*g + f);
q = -(h - w) ./ (w.*g - f);
p = flag_sign * (p + cx);
q = flag_sign * (q - cx);

end

