function p = interpuni(ps,ab)
%INTERPUNI   Piecewise-uniform interpolation.
%
%   P = INTERPUNI(PS,[A B]) constructs a piecewise-polynomial interpolant P
%   whose graph passes through given points on a piecewise-uniform grid.
%
%   Inputs:
%
%     PS: vertical coordinates of interpolation points
%
%     AB = [A B]: interval endpoints
%
%   Outputs:
%
%     P: the interpolating piecewise-polynomial function
%
%   Suppose PS is (M+1)-by-N. Let XS be the M-by-N piecewise-uniform grid
%   on [A,B]. P is constructed so that P(XS(I,J)) = PS(I,J) and so that P
%   is a polynomial of degree at most M on each subinterval in the grid.
%
%   Example:
%
%     f = @(x) sin(x); a = -pi; b = pi;
%     m = 2; n = 4;
%     ps = sampleuni(f,[a b],m,n);
%     p = interpuni(ps,[a b]);
%     newfig;
%     plotfun(f,[a b]);
%     plotpartition([a b],n);
%     plotsample(griduni([a b],m,n),ps);
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('interpuni',ps,ab);
% order subintervals from left to right on the real line
if ab(1)>ab(2)
  ps = fliplr(flipud(ps));
  ab = ab(end:-1:1);
end
% construct piecewise-polynomial interpolant
a = ab(1);
b = ab(2);
m = size(ps,1)-1;
n = size(ps,2);
l = (b-a)/n;
[xs,ws] = griduni([0 1],m,1);
p = @(x) interpuni_(xs,ws,ps,a,b,l,n,x);

