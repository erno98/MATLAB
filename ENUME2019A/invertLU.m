function inverted = invertLU(A)
    
    s = size(A);
    s = s(1);
    
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

