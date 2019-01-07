function [xs,ws] = gridchebpw(asbs,m)
%GRIDCHEBPW   Piecewise-Chebyshev grid.
%
%   [XS,WS] = GRIDCHEBPW(ASBS,M) constructs a piecewise-Chebyshev grid. The
%   subintervals do not have to be of equal width.
%
%   Inputs:
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     [AJ,BJ] in the piecewise-defined grid
%
%     M: the degree of the grid on each subinterval
%
%   Outputs:
%
%     XS: an (M+1)-by-N matrix whose Jth column contains the nodes of a
%     Chebyshev grid of degree M on the Jth subinterval
%
%     WS: an (M+1)-by-1 vector containing the barycentric weights for a
%     Chebyshev grid of degree M
%
%   Example:
%
%     asbs = [ 0 1 3 6; 1 3 6 10 ]
%     m = 6;
%     [xs,ws] = gridchebpw(asbs,m)
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('gridchebpw',asbs,m);
n = size(asbs,2);
xs = nan(m+1,n);
for j = 1:n
  xs(:,j) = gridcheb(asbs(:,j),m);
end
[~,ws] = gridcheb([-1 1],m);

