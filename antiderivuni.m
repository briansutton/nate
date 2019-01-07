function ps = antiderivuni(pa,qs,ab)
%ANTIDERIVUNI   Antiderivative on a piecewise-uniform grid.
%
%   PS = ANTIDERIVUNI(PA,QS,[A B]) computes the antiderivative P of a
%   polynomial Q with initial value P(A) = PA.
%
%   Inputs:
%
%     PA: initial value for the antiderivative
%
%     QS: sample of the integrand Q on an M-by-N piecewise-uniform grid
%     over [A,B]
%
%     AB = [A B]: interval endpoints
%
%   Outputs:
%
%     PS: sample of the antiderivative P on the (M+1)-by-N piecewise-
%     uniform grid over [A,B]
%
%   If Q is a polynomial of degree at most M over each subinterval of the
%   piecewise grid, then the antiderivative is exact, up to roundoff error.
%
%   Example:
%
%     g = @(x) cos(x); a = 0; b = 4*pi;
%     fa = 5;
%     m = 2; n = 6;
%     qs = sampleuni(g,[a b],m,n);
%     ps = antiderivuni(fa,qs,[a b]);
%     p = interpuni(ps,[a b]);
%     f = @(x) 5+sin(x);
%     newfig;
%     plotfun(f,[a b]);
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('antiderivuni',pa,qs,ab);
a = ab(1);
b = ab(2);
m = size(qs,1)-1;
n = size(qs,2);
K = (b-a)/n*intmatrixuni(m);
ps = nan(m+2,n);
pinit = pa;
for j = 1:n
  ps(:,j) = pinit*ones(m+2,1)+K*qs(:,j);
  pinit = ps(end,j);
end

