function tt2j = eds_get_pole(aa, aa1, ha, j)

    ttj = mod(j / sqrt(2), 1);

    % Starting point, assuming a logarithmic distribution is a good guess
    p = polyfit([ log(ha), log(sqrt(ha)) ], [ 0, aa(sqrt(ha)) ], 1);    
    x = exp((ttj - p(2)) ./ p(1)); x = max(ha, min(x, 1));
    
    try
        tt2j = newton(@(y) aa(y) - ttj, aa1, x);
    catch
        tt2j = nan;
    end
    
    if isnan(tt2j) || ~isreal(tt2j)
        tt2j = bisec(ha, 1, -ttj, @(y) aa(y) - ttj);
    end
    
    tt2j = real(sqrt(tt2j));
end