function p = interpcheb(ps,ab)
%INTERPCHEB   Chebyshev interpolation.
%
%   P = INTERPCHEB(PS,[A B]) constructs a polynomial interpolant P with
%   given values on a Chebyshev grid.
%
%   Inputs:
%
%     PS: vertical coordinates of the interpolation points
%
%     AB = [A B]: interval endpoints
%
%   Outputs:
%
%     P: the interpolating polynomial
%
%   Suppose PS is (M+1)-by-1, and let XS be the Chebyshev grid of degree M
%   on [A,B]. Then P is the unique polynomial of degree at most M whose
%   graph intersects the points (XS,PS).
%
%   Example:
%
%     f = @(x) sin(x); a = -pi; b = pi;
%     m = 6;
%     ps = samplecheb(f,[a b],m);
%     p = interpcheb(ps,[a b]);
%     newfig;
%     plotfun(f,[a b]);
%     plotsample(gridcheb([a b],m),ps);
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('interpcheb',ps,ab);
m = size(ps,1)-1;
[xs,ws] = gridcheb(ab,m);
p = @(x) interp_(xs,ws,ps,x);

