function I = definitecheb(qs,ab)
%DEFINITECHEB   Definite integral on a Chebyshev grid.
%
%   I = DEFINITECHEB(QS,[A B]) computes the definite integral of a
%   polynomial Q over [A,B].
%
%   Inputs:
%
%     QS: sample of the integrand on a Chebyshev grid of degree M over
%     [A,B]
%
%     AB = [A B]: limits of integration
%
%   Outputs:
%
%     I: the definite integral
%
%   If Q is a polynomial of degree at most M, then the result is exact, up
%   to roundoff error.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('definitecheb',qs,ab);
a = ab(1);
b = ab(2);
m = size(qs,1)-1;
K = intmatrixcheb(m);
lambdas = K(end,:)';
I = (b-a)/2*lambdas'*qs;

