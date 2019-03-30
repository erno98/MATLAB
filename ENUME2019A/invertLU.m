function inverted = invertLU(A)
    
    s = size(A);
    s = s(1);

    % LU factorization
    [L, U] = lu(A);
    
    % creating identity matrix
    I = eye(s);
    
    y = zeros(s);
    
    for i = 1:s
         y(:, i) = L*U \ I(:, i);
    end

    inverted = y;
   
end

