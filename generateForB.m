function [outputMatrix, outputVector] = generateForB(size)

%for the clarification of output format is set to short
format short

% generate matrix and vector of given sizes
% and filling them up with zeroes
outputMatrix = zeros(size, size);
outputVector = zeros(size, 1);


% filling up matrix with given parameters from formula
for i = 1:size
   for j = 1:size
       outputMatrix(i,j) = 5/(8*(i+j+1));   
   end
end

%filling up vector
for i = 1:size
   if(mod(i,2) ~= 0)
       outputVector(i) = 9/(2*i);
   end
end
