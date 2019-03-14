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

