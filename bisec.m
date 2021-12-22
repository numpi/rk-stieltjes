function x = bisec(a, b, fa, f)
%

tol = eps;

while b - a > tol
    xm = (a+b)/2;
    fm = f(xm);
    
    if fm == 0
        x = xm;
        return;
    end
    
    if fa * fm < 0
        b = xm;
    else
        a = xm;
        fa = fm;
    end
end

x = (b+a)/2;

