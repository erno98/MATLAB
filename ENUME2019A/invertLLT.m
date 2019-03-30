function inverted = invertLLT(A)
    
    s = size(A);
    s = s(1);

    % LLT factorization
    L = chol(A);
    
    % creating identity matrix
    I = eye(s);
    
    y = zeros(s);
    
    for i = 1:s
         y(:, i) = L'*L \ I(:, i);
    end
    
   	inverted = y;
    
    
end

