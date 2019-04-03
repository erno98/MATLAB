clear all;
close all;

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
    
end

% RMSE on x
figure(1)
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
figure(2)
loglog(xs, ME3(1, :), 'r', xs, ME3(2, :), '--r')
hold on
loglog(xs, ME10(1, :),'b',  xs, ME10(2, :), '--b')
loglog(xs, norm10_e(1, :), 'k', xs, norm10_e(2, :), '--k')
loglog(xs, ME20(1, :), 'm', xs, ME20(2, :), '--m')
title('ME on x')
xlabel('x')
ylabel('ME')
legend('LU N=3', 'LLT N=3', 'LU N=10', 'LLT N=10', 'LU N=20', 'LLT N=20', 'Location', 'northwest')
hold off

% comparison with norm
figure(3)
diff_e_LU = abs((norm3_e(1, :) - ME3(1, :)) + (norm10_e(1, :) - ME10(1, :)) + (norm20_e(1, :) - ME20(1, :))) / 3;
diff_e_LLT = abs((norm3_e(2, :) - ME3(2, :)) + (norm10_e(2, :) - ME10(2, :)) + (norm20_e(2, :) - ME20(2, :))) / 3;
diff_inf_LU = abs((norm3_inf(1, :) - RMSE3(1, :)) + (norm10_inf(1, :) - RMSE10(1, :)) + (norm20_inf(1, :) - RMSE20(1, :))) / 3;
diff_inf_LLT = abs((norm3_inf(2, :) - RMSE3(2, :)) + (norm10_inf(2, :) - RMSE10(2, :)) + (norm20_inf(2, :) - RMSE20(2, :))) / 3;
loglog(xs, diff_e_LU, 'ko')
hold on
loglog(xs, diff_e_LLT, 'ro')
loglog(xs, diff_inf_LU, 'k.')
loglog(xs, diff_inf_LLT, 'r.')
xlabel('x')
ylabel('mean difference')
legend('norm-ME LU', 'norm-ME LLT', 'norm.inf-RMSE LU', 'norm.inf-RMSE LLT')
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