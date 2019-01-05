function [J,K] = indefinitecheb(m)
%INDEFINITECHEB   Indefinite integral on Chebyshev grids.
%
%   [J,K] = INDEFINITECHEB(M) computes a pair of matrices J, K for testing
%   if one polynomial is an antiderivative of another.
%
%   Inputs:
%
%     M: degree of a grid
%
%   Outputs:
%
%     J: matrix for evaluating P(X) - P(A) on a Chebyshev grid
%
%     K: matrix for evaluating the integral of Q(X) on a Chebyshev grid
%
%   If PS is a sample of a degree-(M+1) polynomial P on a Chebyshev grid
%   over [A,B] and QS is a sample of a degree-M polynomial Q on a Chebyshev
%   grid over the same interval, then P is an antiderivative of Q if and
%   only if J*PS = (B-A)/2*K*QS.
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('indefinitecheb',m);
J = [ -ones(m+1,1) eye(m+1) ];
K = intmatrixcheb(m);
K = K(2:end,:);

