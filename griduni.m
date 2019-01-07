function [xs,ws] = griduni(ab,m,n)
%GRIDUNI   Piecewise-uniform grid.
%
%   [XS,WS] = GRIDUNI([A B],M,N) constructs an M-by-N piecewise-uniform
%   grid XS on [A,B] and the corresponding barycentric weights WS.
%
%   Inputs:
%
%     AB = [A B]: interval endpoints
%
%     M: the degree of the grid on each subinterval
%
%     N: the number of subintervals (default: N = 1)
%
%   Outputs:
%
%     XS: an (M+1)-by-N matrix containing the nodes of the grid
%
%     WS: an (M+1)-by-1 vector containing the barycentric weights for a
%     uniform grid
%
%   Example:
%
%     a = 0; b = 2;
%     m = 2; n = 4;
%     [xs,ws] = griduni([a b],m,n)
%
%   Copyright 2019 Brian Sutton

narginchk(2,3);
if nargin<3, n = 1; end
natecheck('griduni',ab,m,n);
a = ab(1);
b = ab(2);
if m==0
  % compute nodes and barycentric weights for degree 0
  l = (b-a)/n;
  xs = linspace(a+l/2,b-l/2,n);
  ws = [ 1 ];
else
  % compute nodes for degree >= 1
  xs = linspace(a,b,m*n+1);
  xs = [ reshape(xs(1:end-1),m,n); xs(m+1:m:end) ];
  % compute barycentric weights for degree >= 1
  ws = nan(m+1,1);
  ws(1) = 1;
  for i = 1:floor(m/2)
    ws(i+1) = -ws(i)*(m-i+1)/i;
  end
  k = ceil(m/2);
  ws(end:-1:end-k+1) = (-1)^m*ws(1:k);
end

