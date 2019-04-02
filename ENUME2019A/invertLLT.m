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

