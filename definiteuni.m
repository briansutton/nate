function I = definiteuni(qs,ab)
%DEFINITEUNI   Definite integral on a piecewise-uniform grid.
%
%   I = DEFINITEUNI(QS,[A B]) computes the definite integral of a
%   piecewise-polynomial Q over [A,B].
%
%   Inputs:
%
%     QS: sample of the integrand Q on an M-by-N piecewise-uniform grid
%     over [A,B]
%
%     AB = [A B]: limits of integration
%
%   Outputs:
%
%     I: the definite integral
%
%   The result is exact if Q is a polynomial of degree at most M on each
%   subinterval of the piecewise-uniform grid.
%
%   Example:
%
%     g = @(x) exp(x/2); a = -1; b = 1;
%     m = 2; n = 100;
%     qs = sampleuni(g,[a b],m,n);
%     definiteuni(qs,[a b])
%     2*exp(b/2)-2*exp(a/2)
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('definiteuni',qs,ab);
a = ab(1);
b = ab(2);
m = size(qs,1)-1;
n = size(qs,2);
K = intmatrixuni(m);
lambdas = K(end,:)';
I = (b-a)/n*sum(lambdas'*qs);

