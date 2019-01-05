function I = definitegen(qs,ab,xs,xs_)
%DEFINITEGEN   Definite integral on a given grid.
%
%   DEFINITEGEN computes the definite integral of a polynomial Q over
%   [A,B].
%
%   Inputs:
%
%     QS: sample of the integrand Q on the grid XS_
%
%     AB = [A B]: limits of integration
%
%     XS: any grid of degree M+1
%
%     XS_: grid on which the integrand Q is sampled
%
%   Outputs:
%
%     I: definite integral of Q over [A,B]
%
%   In exact arithmetic, the result is exact if Q is a polynomial of degree
%   at most M.
%
%   Copyright 2019 Brian Sutton

narginchk(4,4);
natecheck('definitegen',qs,ab,xs,xs_);
a = ab(1);
b = ab(2);
K = intmatrixgen(xs,xs_,a);
Eb = evalmatrixgen(b,xs);
lambdas = (Eb*K)';
I = lambdas'*qs;

