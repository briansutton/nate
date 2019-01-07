function x = gepp(A,b)
%GEPP   Solve a linear system using Gaussian elimination with partial
%pivoting.
%
%   X = GEPP(A,B) solves
%
%     A*X = B
%
%   for X, in which A is a square matrix and B is a column vector.
%
%   Inputs:
%
%     A: an N-by-N matrix
%
%     B: an N-by-1 vector
%
%   Outputs:
%
%     X: the solution to A*X = B
%
%   The computed solution may be inaccurate if A is ill conditioned or
%   meaningless if A*X = B does not have a solution.
%
%   Example:
%
%     A = [ 1 2; 3 4 ]
%     b = [4; 10]
%     x = gepp(A,b)
%     % check
%     A*x
%     b
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('gepp',A,b);
n = size(A,1);
U = A;
y = b;
% for each column,
for j = 1:n-1
  % find pivot
  [~,k] = max(abs(U(j:n,j)));
  k = j+k-1;
  if k>j
    % swap rows
    U([j k],j:n) = U([k j],j:n);
    y([j k]) = y([k j]);
  end
  % eliminate entries below pivot
  l = U(j+1:n,j)/U(j,j);
  U(j+1:n,j) = 0;
  U(j+1:n,j+1:n) = U(j+1:n,j+1:n)-l*U(j,j+1:n);
  % apply same operation to right-hand side
  y(j+1:n) = y(j+1:n)-l*y(j);
end
% solve upper-triangular system
x = backsubstitute(U,y);

