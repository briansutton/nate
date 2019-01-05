function [J,K] = indefiniteuni(m)
%INDEFINITEUNI   Indefinite integral on uniform grids.
%
%   [J,K] = INDEFINITEUNI(M) computes a pair of matrices J, K for testing
%   if one polynomial is an antiderivative of another.
%
%   Inputs:
%
%     M: degree of a grid
%
%   Outputs:
%
%     J: matrix for evaluating P(X) - P(A) on a uniform grid
%
%     K: matrix for evaluating the integral of Q(X) on a uniform grid
%
%   If PS is a sample of a degree-(M+1) polynomial P on a uniform grid over
%   [A,B] and QS is a sample of a degree-M polynomial Q on a uniform grid
%   over the same interval, then P is an antiderivative of Q if and only if
%   J*PS = (B-A)*K*QS.
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('indefiniteuni',m);
J = [ -ones(m+1,1) eye(m+1) ];
K = intmatrixuni(m);
K = K(2:end,:);

