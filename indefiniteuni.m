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
%   Example:
%
%     m = 2;
%     [J,K] = indefiniteuni(m);
%     q = @(x) 5*x^2+2*x-3;
%     p = @(x) (5/3)*x^3+x^2-3*x+4;
%     a = -2; b = 2;
%     qs = sampleuni(q,[a b],m);
%     ps = sampleuni(p,[a b],m+1);
%     J*ps
%     (b-a)*K*qs
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('indefiniteuni',m);
J = [ -ones(m+1,1) eye(m+1) ];
K = intmatrixuni(m);
K = K(2:end,:);

