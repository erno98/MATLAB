close all
clear all

% for determinant = 0, x has to be 0, thus:
fun = @f;
% alpha has to be positive, starting from 0
alpha = abs(fzero(fun, 0));
x = sqrt(alpha^2 +1/2) -1;

% matrix generation
A1 = generate_matrix(3, x);
A2 = generate_matrix(10, x);
A3 = generate_matrix(20, x);

% dependances
num = 700;
alphas = linspace(alpha - 0.01, alpha + 0.01, num);
dets3 = zeros(num, 1);
conds3 = zeros(num, 1);

dets10 = zeros(num, 1);
conds10 = zeros(num, 1);

dets20 = zeros(num, 1);
conds20 = zeros(num, 1);



for i = 1:num
    A3 = generate_matrix(3, sqrt(alphas(i)^2 +1/2) -1);
    A10 = generate_matrix(10, sqrt(alphas(i)^2 +1/2) -1);
    A20 = generate_matrix(20, sqrt(alphas(i)^2 +1/2) -1);
    
    dets3(i) =  det(A3);
    conds3(i) = cond(A3);
    
    dets10(i) =  det(A10);
    conds10(i) = cond(A10);
    
    dets20(i) =  det(A20);
    conds20(i) = cond(A20);
end

% det(A) on alpha
figure(1)
semilogy(alphas, dets3, alphas, dets10, alphas, dets20)
legend("N = 3", "N = 10", "N = 20")
title("det(A) on alpha")
xlabel("alpha")
ylabel("det(A)")

% cond(A) on alpha
figure(2)
semilogy(alphas, conds3, alphas, conds10, alphas, conds20)
legend("N = 3", "N = 10", "N = 20")
title("cond(A) on alpha")
xlabel("alpha")
ylabel("cond(A)")

function x = f(alpha)
    x = sqrt(alpha^2 + 1/2) - 1;
end
