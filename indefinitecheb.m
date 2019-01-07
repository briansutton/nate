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
%   Example:
%
%     m = 2;
%     [J,K] = indefinitecheb(m);
%     q = @(x) 5*x^2+2*x-3;
%     p = @(x) (5/3)*x^3+x^2-3*x+4;
%     a = -2; b = 2;
%     qs = samplecheb(q,[a b],m);
%     ps = samplecheb(p,[a b],m+1);
%     J*ps
%     ((b-a)/2)*K*qs
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('indefinitecheb',m);
J = [ -ones(m+1,1) eye(m+1) ];
K = intmatrixcheb(m);
K = K(2:end,:);

