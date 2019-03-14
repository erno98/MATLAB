% setting up initial parameters
precision = 10e-7;
step = 10e-7;
[alpha_1, alpha_set_1, det_set_1, cond_set_1] = find_alpha(3, precision, step);
[alpha_10, alpha_set_10, det_set_10, cond_set_10] = find_alpha(10, precision, step);
[alpha_20, alpha_set_20, det_set_20, cond_set_20] = find_alpha(20, precision, step);

% -------------------------------------------

figure('Name', 'Det(A) on alpha')
subplot(1,3,1);
plot(alpha_set_1, det_set_1)
title('N = 3');
xlabel('alpha')
ylabel('det(A)')

subplot(1,3,2);
plot(alpha_set_10, det_set_10)
title('N = 10');
xlabel('alpha')
ylabel('det(A)')

subplot(1,3,3);
plot(alpha_set_20, det_set_20)
title('N = 20');
xlabel('alpha')
ylabel('det(A)')

% -------------------------------------------

figure('Name', 'cond(A) on alpha') 
subplot(1,3,1);
plot(alpha_set_1, cond_set_1)
title('N = 3');
xlabel('alpha')
ylabel('cond(A)')

subplot(1,3,2);
plot(alpha_set_10, cond_set_10)
title('N = 10');
xlabel('alpha')
ylabel('cond(A)')

subplot(1,3,3);
plot(alpha_set_20, cond_set_20)
title('N = 20');
xlabel('alpha')
ylabel('cond(A)')


