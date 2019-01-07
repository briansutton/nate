function [L,U,P] = lupp(A)
%LUPP   LU factorization with partial pivoting.
%
%   [L,U,P] = LUPP(A) computes an LU factorization
%
%     P*A = L*U,
%
%   in which P is a permutation matrix.
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
%     P: permutation matrix
%
%   Example:
%
%     A = [ 1 2; 3 4 ]
%     [L,U,P] = lupp(A)
%     % check
%     P*A-L*U
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('lupp',A);
n = size(A,1);
sigma = 1:n;
L = eye(n);
U = A;
% for each column,
for j = 1:n-1
  % find pivot
  [~,k] = max(abs(U(j:n,j)));
  k = j+k-1;
  if k>j
    % swap rows
    sigma([j k]) = sigma([k j]);
    L([j k],1:j-1) = L([k j],1:j-1);
    U([j k],j:n) = U([k j],j:n);
  end
  l = U(j+1:n,j)/U(j,j);
  % store and apply elimination operation
  L(j+1:n,j) = l;
  U(j+1:n,j) = 0;
  U(j+1:n,j+1:n) = U(j+1:n,j+1:n)-l*U(j,j+1:n);
end
% construct permutation matrix
P = eye(n);
P = P(sigma,:);

