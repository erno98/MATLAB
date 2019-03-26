function [alpha, alpha_set, det_set, cond_set] = find_alpha(N, precision, step)

% TODO: zrobić alfę f0, przyrównać x do 0 i naprawić ten badziew

    % initial values
    alpha = 0;
    x = sqrt(alpha^2 + 1/2) - 1;
    A = generate_matrix(N, x);
    
    % finding smallest alpha satisfying det(A)=0
    i = 1;
    while det(A) > precision
        alpha = alpha + step;
        x = sqrt(alpha^2 + 1/2) - 1;
        A = generate_matrix(N, x);
        i = i+1;
    end
        
    % outputting alpha, det(A) and cond(A) values for plotting
    i=1;
    a = alpha-0.01;
    while a < alpha+0.01
        alpha_set(i) = a;
        x = sqrt(a^2 + 1/2) - 1;
        A = generate_matrix(N, x);
        det_set(i) = det(A);
        cond_set(i) = cond(A);
        i = i+1;
        a = a + step;
    end

end
