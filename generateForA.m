function [outputMatrix, outputVector] = generateForA(size)

%for the clarification of output format is set to short
format short

% generate matrix and vector of given sizes
% and filling them up with zeroes
outputMatrix = zeros(size, size);
outputVector = zeros(size, 1);


% filling up matrix with 9's and -3's
for i = 1:size
   for j = 1:size
       if( i == j)
            outputMatrix(i, j) = 9;
       end
        if (i == j-1 || i == j+1)
            outputMatrix(i, j) = -3;
        end    
   end
end

%filling up vector
for i = 1:size
   outputVector(i) = (-5 + 0.3*i);
end








