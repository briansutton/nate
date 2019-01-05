function ps = antiderivgen(pa,qs,xs,xs_,a)
%ANTIDERIVGEN   Antiderivative on a given grid.
%
%   PS = ANTIDERIVGEN(PA,QS,XS,XS_,A) computes the antiderivative P of a
%   polynomial Q with initial value P(A) = PA.
%
%   Inputs:
%
%     PA: initial value for the antiderivative
%
%     QS: sample of the integrand Q on the grid XS_
%
%     XS: grid for the antiderivative P (degree M+1)
%
%     XS_: grid for the integrand Q (degree M)
%
%     A: lower limit of integration
%
%   Outputs:
%
%     PS: sample of the antiderivative P on the grid XS
%
%   If Q is a polynomial of degree at most M, then the antiderivative is
%   exact, up to roundoff error.
%
%   Copyright 2019 Brian Sutton

narginchk(5,5);
natecheck('antiderivgen',pa,qs,xs,xs_,a);
m = length(qs)-1;
K = intmatrixgen(xs,xs_,a);
ps = pa*ones(m+2,1)+K*qs;

