k = 0:21;
RMSE3 = zeros(2, length(k));
ME3 = zeros(2, length(k));

RMSE10 = zeros(2, length(k));
ME10 = zeros(2, length(k));

RMSE20 = zeros(2, length(k));
ME20 = zeros(2, length(k));

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
    
    RMSE10(1, i) = norm_2(A10, invertLU(A10));
    RMSE10(2, i) = norm_2(A10, invertLLT(A10));
    ME10(1, i) = norm_inf(A10, invertLU(A10));
    ME10(2, i) = norm_inf(A10, invertLLT(A10));
    
    RMSE20(1, i) = norm_2(A20, invertLU(A20));
    RMSE20(2, i) = norm_2(A20, invertLLT(A20));
    ME20(1, i) = norm_inf(A20, invertLU(A20));
    ME20(2, i) = norm_inf(A20, invertLLT(A20));
    
end

figure(1)
plot(xs, RMSE3(1, :), xs, RMSE3(2, :))
hold on
plot(xs, RMSE10(1, :), xs, RMSE10(2, :))
plot(xs, RMSE20(1, :), xs, RMSE20(2, :))
title('RMSE on x')
xlabel('x')
ylabel('RMSE')
legend('LU N=3', 'LLT N=3', 'LU N=10', 'LLT N=10', 'LU N=20', 'LLT N=20', 'Location', 'northwest')
hold off

disp(RMSE10(1, :) == RMSE10(2,:))

figure(2)
plot(xs, ME3(1, :), xs, ME3(2, :))
hold on
plot(xs, ME10(1, :), xs, ME10(2, :))
plot(xs, ME20(1, :), xs, ME20(2, :))
title('ME on x')
xlabel('x')
ylabel('ME')
legend('LU N=3', 'LLT N=3', 'LU N=10', 'LLT N=10', 'LU N=20', 'LLT N=20', 'Location', 'northwest')
hold off

function rmse = norm_2(A, A_estimate)
    A_to_norm = A*A_estimate - eye(size(A));
    rmse = max(sqrt(max(abs(eig(A_to_norm * A_to_norm.')))));
end

function me = norm_inf(A, A_estimate)
    A_to_norm = A*A_estimate - eye(size(A));
    A_to_norm = abs(A_to_norm);
    me = max(sum(A_to_norm, 2));
end