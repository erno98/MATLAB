
% plotting the function
fplot(@(x) f(x), [-1, 1])

% indicating sequences
s10 = indicate(10);
s20 = indicate(20);
s30 = indicate(30);


y = @(x) (x+1/3)^2 + exp(-x-2);


function sequence = indicate(N)
    sequence = zeros(1, N);
    for n = 1:N
       x = -1 + 2*(n-1)/(N-1);
       sequence(n) = f(x);
    end
end
