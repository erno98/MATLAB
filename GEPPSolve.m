function [Solutions, ResiduumNorm, resCorrected] = GEPPSolve(A,b)

%check whether the matrix is square, and vector is of proper dimension

%creating temporary variables to store inputs size to check them
size1 = size(A);
size2 = size(b);

%checking if the matrix is square
if(size1(1) ~= size1(2))
    error("Given matrix is not square one.")
end

%checking if the vector has proper dimension
if (size2(2) ~= 1) 
    error("Given vector has wrong dimensions.")
end

%checking vector and matrix dimensions to augment them
if(size1(1) ~= size2(1))
   error("Vector and matrix dimension don't match.")
end

%deleting the temporary variables
clear size1 size2;

%-----------------------------------------
%input is correct, proceed to function
%-----------------------------------------

%create augmented matrix with given matrix A and vector b
augmentedMatrix = [A,b];

% 1x1 matrix with size of given square matrix
sizeOfA = size(A,1);

   
for count = 1:sizeOfA
    
%-----------------------------------------
% partial pivoting
%-----------------------------------------
    
   %additional variable not to modify the 'count' from main loop
    i = count;
    for j = count+1 : sizeOfA
        if (augmentedMatrix(j, count) > augmentedMatrix(i, count) )
            %index of row with the biggest value
            i = j;
        end
    end
    %swap rows if necessary
    if( count~=i )
        augmentedMatrix([count i], :) = augmentedMatrix([i count], :);
    end
    
   
%-----------------------------------------
% gaussian elimination
%-----------------------------------------

   %transforming matrix to echelon form
    for j = count+1 : sizeOfA
        L = augmentedMatrix(j, count) / augmentedMatrix(count , count);
        augmentedMatrix(j,:) = augmentedMatrix(j,:) - augmentedMatrix(count,:) * L;
    end
end
    
%-----------------------------------------
% solving the linear equations
%-----------------------------------------


%setting up solutions vector for keeping the answers
Solutions = eye(sizeOfA);

solutionsVector = zeros(1, sizeOfA);

%separating augmented matrix for easier calculations
Aechelon = augmentedMatrix(:, 1:end-1);
bechelon = augmentedMatrix(:, end);

for count = 1: sizeOfA
   k = sizeOfA - count + 1;
   sum = 0;
   for j = k+1 : sizeOfA
      sum = sum + (augmentedMatrix(k,j) * solutionsVector(j)); 
   end
   solutionsVector(k) = (bechelon(k) - sum)/Aechelon(k,k);
end
    
solutionsVector = solutionsVector';
Solutions = [Solutions, solutionsVector];

%-----------------------------------------
% solution error calculation
%-----------------------------------------


%calculating the residuum
residuum = A * solutionsVector - b;

%calculating the Euclidan norm of residuum
norm = 0;
for count = 1 : size(residuum)
    norm = norm + (residuum(count) * residuum(count));
end

%outputting normalized residuum
ResiduumNorm = sqrt(norm);


%residuum correction
deltaX = residuum' * inv(A);
deltaX = deltaX';
x2 = solutionsVector - deltaX;
res2 = A * x2 - b;

%calculating the Euclidan norm of corrected residuum
norm2 = 0;
for count = 1 : size(res2)
    norm2 = norm2 + (res2(count) * res2(count));
end

resCorrected = sqrt(norm2);


end
