function [root, iterations] = NewtonRoot(f, a, b, accuracy)
%function finds one root of non-linear function in given interval

iterations = 0;



%initial value of root
root = (a+b)/2;

%calculating derivative of the function with respect to x
syms x
dfdx = eval(['@(x)' char(diff(f(x)))]);

while(abs(f(root)) > accuracy)
    
    %apply the formula
    root = root - (f(root) / dfdx(root));
    iterations = iterations + 1;

    
end


end

