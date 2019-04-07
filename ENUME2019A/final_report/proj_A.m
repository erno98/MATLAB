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
num = 2000;
alphas = linspace(alpha - 0.01, alpha + 0.01, num);
dets3 = zeros(num, 1);
conds3 = zeros(num, 1);

dets10 = zeros(num, 1);
conds10 = zeros(num, 1);

dets20 = zeros(num, 1);
conds20 = zeros(num, 1);



for i = 1:num
    A3 = generate_matrix(3, f(alphas(i)));
    A10 = generate_matrix(10, f(alphas(i)));
    A20 = generate_matrix(20, f(alphas(i)));
    
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


k = 0:21;
RMSE3 = zeros(2, length(k));
ME3 = zeros(2, length(k));
norm3_e = zeros(2, length(k));
norm3_inf = zeros(2, length(k));

RMSE10 = zeros(2, length(k));
ME10 = zeros(2, length(k));
norm10_e = zeros(2, length(k));
norm10_inf = zeros(2, length(k));

RMSE20 = zeros(2, length(k));
ME20 = zeros(2, length(k));
norm20_e = zeros(2, length(k));
norm20_inf = zeros(2, length(k));

xs = zeros(1, length(k));
means = zeros(6, length(k));

for i = k+1
    xs(i) = 2^(i-1)/300;
    A3 = generate_matrix(3, xs(i));
    A10 = generate_matrix(10, xs(i));
    A20 = generate_matrix(20, xs(i));
    
    RMSE3(1, i) = norm_2(A3, invertLU(A3));
    RMSE3(2, i) = norm_2(A3, invertLLT(A3));
    ME3(1, i) = norm_inf(A3, invertLU(A3));
    ME3(2, i) = norm_inf(A3, invertLLT(A3));
    norm3_e(1, i) = m_norm_euc(A3, invertLU(A3));
    norm3_e(2, i) = m_norm_euc(A3, invertLLT(A3));
    norm3_inf(1, i) = m_norm_inf(A3, invertLU(A3));
    norm3_inf(2, i) = m_norm_inf(A3, invertLLT(A3));
    
    RMSE10(1, i) = norm_2(A10, invertLU(A10));
    RMSE10(2, i) = norm_2(A10, invertLLT(A10));
    ME10(1, i) = norm_inf(A10, invertLU(A10));
    ME10(2, i) = norm_inf(A10, invertLLT(A10));
    norm10_e(1, i) = m_norm_euc(A10, invertLU(A10));
    norm10_e(2, i) = m_norm_euc(A10, invertLLT(A10));
    norm10_inf(1, i) = m_norm_inf(A10, invertLU(A10));
    norm10_inf(2, i) = m_norm_inf(A10, invertLLT(A10));
    
    RMSE20(1, i) = norm_2(A20, invertLU(A20));
    RMSE20(2, i) = norm_2(A20, invertLLT(A20));
    ME20(1, i) = norm_inf(A20, invertLU(A20));
    ME20(2, i) = norm_inf(A20, invertLLT(A20));
    norm20_e(1, i) = m_norm_euc(A20, invertLU(A20));
    norm20_e(2, i) = m_norm_euc(A20, invertLLT(A20));
    norm20_inf(1, i) = m_norm_inf(A20, invertLU(A20));
    norm20_inf(2, i) = m_norm_inf(A20, invertLLT(A20));
    
    means(1, i) = mean(abs(inv(A3) - invertLU(A3)), 'All');
    means(2, i) = mean(abs(inv(A3) - invertLLT(A3)), 'All');
    
    means(3, i) = mean(abs(inv(A10) - invertLU(A10)), 'All');
    means(4, i) = mean(abs(inv(A10) - invertLLT(A10)), 'All');
    
    means(5, i) = mean(abs(inv(A20) - invertLU(A20)), 'All');
    means(6, i) = mean(abs(inv(A20) - invertLLT(A20)), 'All');
    
end



% RMSE on x
figure(3)
loglog(xs, RMSE3(1, :), 'r', xs, RMSE3(2, :), '--r')
% loglog(xs, norm20(1, :), '*k', xs, norm20(2, :), '--k')
hold on
loglog(xs, RMSE10(1, :), 'b', xs, RMSE10(2, :), '--b')
loglog(xs, RMSE20(1, :), 'm', xs, RMSE20(2, :), '--m')
title('RMSE on x')
xlabel('x')
ylabel('RMSE')
legend('LU N=3', 'LLT N=3', 'LU N=10', 'LLT N=10', 'LU N=20', 'LLT N=20', 'Location', 'northwest')
hold off

% ME on x
figure(4)
loglog(xs, ME3(1, :), 'r', xs, ME3(2, :), '--r')
hold on
loglog(xs, ME10(1, :),'b',  xs, ME10(2, :), '--b')
loglog(xs, ME20(1, :), 'm', xs, ME20(2, :), '--m')
title('ME on x')
xlabel('x')
ylabel('ME')
legend('LU N=3', 'LLT N=3', 'LU N=10', 'LLT N=10', 'LU N=20', 'LLT N=20', 'Location', 'northwest')
hold off

% check norm 
norm_euc_diff_3_LU = mean(abs(ME3(1, :) - norm3_e(1, :)));
norm_euc_diff_10_LU = mean(abs(ME10(1, :) - norm10_e(1, :)));
norm_euc_diff_20_LU = mean(abs(ME20(1, :) - norm20_e(1, :)));

norm_inf_diff_3_LU = mean(abs(RMSE3(1, :) - norm3_inf(1, :)));
norm_inf_diff_10_LU = mean(abs(RMSE10(1, :) - norm10_inf(1, :)));
norm_inf_diff_20_LU = mean(abs(RMSE20(1, :) - norm20_inf(1, :)));

norm_euc_diff_3_LLT = mean(abs(ME3(2, :) - norm3_e(2, :)));
norm_euc_diff_10_LLT = mean(abs(ME10(2, :) - norm10_e(2, :)));
norm_euc_diff_20_LLT = mean(abs(ME20(2, :) - norm20_e(2, :)));

norm_inf_diff_3_LLT = mean(abs(RMSE3(2, :) - norm3_inf(2, :)));
norm_inf_diff_10_LLT = mean(abs(RMSE10(2, :) - norm10_inf(2, :)));
norm_inf_diff_20_LLT = mean(abs(RMSE20(2, :) - norm20_inf(2, :)));

c = categorical({'euclidan-LU', 'infinite-LU', 'euclidan-LLT', 'infinite-LLT'});
norms = [norm_euc_diff_3_LU, norm_euc_diff_10_LU, norm_euc_diff_20_LU; ...
    norm_inf_diff_3_LU, norm_inf_diff_10_LU, norm_inf_diff_20_LU; ...
    norm_euc_diff_3_LLT, norm_euc_diff_10_LLT, norm_euc_diff_20_LLT; ...
    norm_inf_diff_3_LLT, norm_inf_diff_10_LLT, norm_inf_diff_10_LLT
    ];

figure(5)
bar(c, norms)
hold on
title("Average difference between matlab's norms and created ones")
ylabel("norm difference")
xlabel('norm-decomposition')
legend('N=3', 'N=10', 'N=20', 'Location', 'northwest')
grid on
hold off

% comparison to inv
figure(6)
means_to_bar = [mean(means(1, :)), mean(means(2, :)); ...
                mean(means(3, :)), mean(means(4, :)); ...
                mean(means(5, :)), mean(means(6, :))
];

c = categorical({'N=3', 'N=10', 'N=20'});
c = reordercats(c, {'N=3', 'N=10', 'N=20'});

bar(c, means_to_bar)
hold on
set(gca,'YScale','log')
grid on
title("average difference between inverse and inv()")
ylabel('difference')
legend('LU', 'LLT', 'Location', 'northeast')
hold off

function rmse = norm_2(A, A_estimate)
    A_to_norm = A*A_estimate - eye(size(A));
    rmse = max(sqrt(max(abs(eig(A_to_norm * A_to_norm.')))));
end

function m = m_norm_euc(A, A_estimate)
    A_to_norm = A*A_estimate - eye(size(A));
    m = norm(A_to_norm);
end

function m = m_norm_inf(A, A_estimate)
    A_to_norm = A*A_estimate - eye(size(A));
    m = norm(A_to_norm, Inf);
end

function me = norm_inf(A, A_estimate)
    A_to_norm = A*A_estimate - eye(size(A));
    A_to_norm = abs(A_to_norm);
    me = max(sum(A_to_norm, 2));
end

function matrix = generate_matrix(N, x)

    % creating matrix of given size filled with 0s
    matrix = zeros(N, N);
    
    % first row and column of the matrix
    matrix(1,1) = x^2;
    for i = 2:N
       matrix(1, i) = (3*x)/(( (-1)^(i-1) )* 2); 
       matrix(i, 1) = (3*x)/(( (-1)^(i-1) )* 2); 
    end
    
    % other elements of the matrix
    for i = 2:N
       for j = 2:N
            matrix(i, j) = (9*j)/(4*(-1)^(j-i));
            matrix(j, i) = matrix(i, j);
       end
    end
end

function x = f(alpha)
    x = sqrt(alpha^2 + 1/2) - 1;
end

function inverted = invertLLT(A)
    
    s = size(A);
    s = s(1);

    % LLT factorization
    LT = chol(A);
    L = transpose(LT);
    
    % creating identity matrix
    I = eye(s);
    
    inverted = zeros(s);
    X = zeros(s);
    
    for k=1:s
        X(1,k) = I(1,k)./L(1,1);
        for m=2:s
            X(m,k) = (I(m,k)-L(m,1:m-1)*X(1:m-1,k))./L(m,m);
        end
        
        inverted(s,k) = X(s,k)./LT(s,s);
        
        for m=s-1:-1:1
            inverted(m,k) = (X(m,k)-LT(m,m+1:s)*inverted(m+1:s,k))./LT(m,m);
        end
    end
    
end

function inverted = invertLU(A)
    
    s = size(A);
    s = s(1);
    
    % LU factorization
    [L, U, P] = lu(A);
   
    % creating identity matrix
    I = eye(s);
    
    inverted = zeros(s);
    X = zeros(s);

    
    for k=1:s
        X(1,k) = I(1,k);
        for m=2:s
            X(m,k) = (I(m,k)-L(m,1:m-1)*X(1:m-1,k));
        end
        
        inverted(s,k) = X(s,k)./U(s,s);
        
        for m=s-1:-1:1
            inverted(m,k) = (X(m,k)-U(m,m+1:s)*inverted(m+1:s,k))./U(m,m);
        end
    end
   
    inverted = inverted*P;

    
end

