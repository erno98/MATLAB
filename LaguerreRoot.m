function [root, iterations] = LaguerreRoot(P, xk, accuracy)


% First derivative
dP = polyder(P);

% Second derivative
dP2 = polyder(dP);

% Initial values

n = length(P) - 1;
iterations = 0;

while(polyval(P, xk) > accuracy)
    
    %squre root value in the denominator
    denomSqrt = sqrt( (n-1) * ( (n-1) * (polyval(dP, xk)^2) - n*polyval(P, xk) *polyval(dP2, xk)));
    
    %choosing largest absolute value for the denominator 
    if(abs(polyval(dP, xk) + denomSqrt) > abs(polyval(dP, xk) - denomSqrt))
        denom = polyval(dP, xk) + denomSqrt;
    else 
        denom = polyval(dP, xk) - denomSqrt;
    end
    

    
    %applying the formula
    xk = xk - (n*polyval(P,xk) / denom);
    iterations = iterations + 1;
end

root = xk;

end

