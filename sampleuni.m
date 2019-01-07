function ys = sampleuni(f,ab,m,n)
%SAMPLE   Sample on an M-by-N piecewise-uniform grid.
%
%   YS = SAMPLE(F,[A B],M,N) evaluates a function at each node in a
%   piecewise-uniform grid.
%
%   Inputs:
%
%     F: function
%
%     AB = [A B]: endpoints of an interval
%
%     M: degree
%
%     N: number of subintervals (default: N = 1)
%
%   Outputs:
%
%     YS: matrix of function values
%
%   An M-by-N piecewise-uniform grid is constructed on [A,B], and then the
%   function F is evaluated at each node to produce YS.
%
%   Example:
%
%     f = @(x) cos(x); a = 0; b = 2*pi;
%     m = 2; n = 4;
%     xs = griduni([a b],m,n)
%     ps = sampleuni(f,[a b],m,n)
%     newfig;
%     plotfun(f,[a b]);
%     plotpartition([a b],n);
%     plotsample(xs,ps);
%
%   Copyright 2019 Brian Sutton

narginchk(3,4);
if nargin<4, n = 1; end
natecheck('sampleuni',f,ab,m,n);
a = ab(1);
b = ab(2);
if m==0
  l = (b-a)/n;
  ys = arrayfun(f,linspace(a+l/2,b-l/2,n));
else
  ys = arrayfun(f,linspace(a,b,m*n+1));
  ys = [ reshape(ys(1:end-1),m,n); ys(m+1:m:end) ];
end


