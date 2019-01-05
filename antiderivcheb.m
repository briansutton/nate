function ps = antiderivcheb(pa,qs,ab)
%ANTIDERIVCHEB   Antiderivative on a Chebyshev grid.
%
%   PS = ANTIDERIVCHEB(PA,QS,[A B]) computes the antiderivative P of a
%   polynomial Q with initial value P(A) = PA.
%
%   Inputs:
%
%     PA: initial value for the antiderivative
%
%     QS: sample of the integrand Q on a Chebyshev grid of degree M over
%     [A,B]
%
%     AB = [A B]: interval endpoints
%
%   Outputs:
%
%     PS: sample of the antiderivative on the Chebyshev grid of degree M+1
%     over [A,B]
%
%   If Q is a polynomial of degree at most M, then the antiderivative is
%   exact, up to roundoff error.
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('antiderivcheb',pa,qs,ab);
a = ab(1);
b = ab(2);
m = length(qs)-1;
K = (b-a)/2*intmatrixcheb(m);
ps = pa*ones(m+2,1)+K*qs;


