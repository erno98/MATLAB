function [root, iterations] = BisectionRoot(f, a, b, accuracy)
%function finds one root of non-linear function in given interval

%initial value of root
div_point = (a + b) / 2;
iterations = 0;


while(abs(f(div_point)) > accuracy)

    
    root = div_point;
    
    %choosing new interval
    if(f(a)*f(div_point) < 0)
        b = div_point;
    elseif (f(div_point)*f(b) < 0)
        a = div_point;
    end
    
    %dividing the interval into two equal ones
     div_point = (a + b) / 2;
     
     iterations = iterations + 1;
    
end


end

