function z = zerosuni(ps,ab,tol)
%ZEROSUNI   Zeros of a piecewise-polynomial interpolant on a piecewise-
%uniform grid.
%
%   Z = ZEROSUNI(PS,[A B],TOL) computes the zeros of an interpolant on a
%   piecewise-uniform grid.
%
%   Inputs:
%
%     PS: sample of a piecewise-polynomial function on a piecewise-uniform
%     grid over [A,B]
%
%     AB = [A B]: endpoints of the interval
%
%     TOL: tolerance used in testing if a computed zero is close enough to
%     the given interval (default: TOL = 1e-6)
%
%   Outputs:
%
%     Z: zeros of the interpolant in the interval [A,B]
%
%   If the grid is M-by-N and the polynomial pieces are of degree at most
%   M, then the zeros are exact (in theory).
%
%   Example:
%
%     f = @(x) cos(x); a = 0; b = 10;
%     m = 3; n = 100;
%     ps = sampleuni(f,[a b],m,n);
%     z = zerosuni(ps,[a b])
%     newfig;
%     plotfun(f,[a b]);
%     plot(z,zeros(size(z)),'*');
%
%   Copyright 2019 Brian Sutton

narginchk(2,3);
if nargin<3, tol = 1e-6; end
natecheck('zerosuni',ps,ab,tol);
m = size(ps,1)-1;
n = size(ps,2);
a = ab(1);
b = ab(2);
l = (b-a)/n;
[xs,ws] = griduni([0 1],m);
z = [];
for j = 1:n
  % find zeros in subinterval
  e = zeros_(xs,ws,ps(:,j));
  e = e(abs(imag(e))<=tol&real(e)>=-tol&real(e)<=1+tol);
  % append zeros to running list
  z = [ z; a+l*(j-1+e) ];
end
z = sort(z);

