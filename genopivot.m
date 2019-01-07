function x = genopivot(A,b)
%GENOPIVOT   Solve a linear system using Gaussian elimination without
%pivoting.
%
%   X = GENOPIVOT(A,B) solves
%
%     A*X = B
%
%   for X, in which A is a square matrix and B is a column vector.
%
%   Inputs:
%
%     A: an N-by-N matrix
%
%     B: an N-by-1 column vector
%
%   Outputs:
%
%     X: the solution to A*X = B
%
%   This routine implements Gaussian elimination without pivoting. It may
%   produce garbage if any pivot equals zero.
%
%   Example:
%
%     A = [ 4 3; 2 1 ]
%     b = [10; 4]
%     x = genopivot(A,b)
%     % check
%     A*x
%     b
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('genopivot',A,b);
n = size(A,1);
U = A;
y = b;
% for each column,
for j = 1:n-1
  % eliminate entries below pivot
  v = U(j+1:end,j)/U(j,j);
  U(j+1:end,j) = 0;
  U(j+1:end,j+1:end) = U(j+1:end,j+1:end)-v*U(j,j+1:end);
  % apply same operation to right-hand side
  y(j+1:end) = y(j+1:end)-y(j)*v;
end
% solve upper-triangular system
x = backsubstitute(U,y);

