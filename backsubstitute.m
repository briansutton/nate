function x = backsubstitute(U,b)
%BACKSUBSTITUTE  Solve an upper-triangular linear system by back
%substitution.
%
%   X = BACKSUBSTITUTE(U,B) computes the solution X to
%
%     U*X = B,
%
%   in which U is a square upper-triangular matrix and B is a column
%   vector.
%
%   Inputs:
%
%     U: an N-by-N upper-triangular matrix
%
%     B: an N-by-1 column vector
%
%   Outputs:
%
%     X: the solution to U*X = B
%
%   Example:
%
%     U = [ 2 3; 0 5 ]
%     b = [1; 5]
%     x = backsubstitute(U,b)
%     % check
%     U*x
%     b
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('backsubstitute',U,b);
n = size(U,1);
x = nan(n,1);
for i = n:-1:1
  x(i) = (b(i)-U(i,i+1:n)*x(i+1:n))/U(i,i);
end

