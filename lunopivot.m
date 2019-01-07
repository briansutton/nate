function [L,U] = lunopivot(A)
%LUNOPIVOT   LU factorization without pivoting.
%
%   [L,U] = LUNOPIVOT(A) computes an LU factorization
%
%     A = L*U.
%
%   The routine fails if any pivot equals zero.
%
%   Inputs:
%
%     A: a square matrix
%
%   Outputs:
%
%     L: unit lower-triangular factor
%
%     U: upper-triangular factor
%
%   Example:
%
%     A = [ 4 3; 2 1 ]
%     [L,U] = lunopivot(A)
%     % check
%     A-L*U
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('lunopivot',A);
n = size(A,1);
L = eye(n);
U = A;
% for each column,
for j = 1:n-1
  v = U(j+1:n,j)/U(j,j);
  % store elimination operation
  L(j+1:n,j) = v;
  % apply elimination operation
  U(j+1:n,j) = 0;
  U(j+1:n,j+1:n) = U(j+1:n,j+1:n)-v*U(j,j+1:n);
end

