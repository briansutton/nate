function D = diffmatrixsquare_(xs,ws)
%DIFFMATRIXSQUARE_   Square differentiation matrix on a single grid.
%
%   D = DIFFMATRIXSQUARE_(XS,WS) computes a square differentiation matrix
%   using identical source and destination grids.
%
%   Inputs:
%
%     XS: a grid of degree M
%
%     WS: barycentric weights for the grid XS
%
%   Outputs:
%
%     D: the differentiation matrix
%
%   If PS is a sample of a degree-M polynomial P on the grid XS, then QS =
%   D*PS is a sample of Q = P' on the same grid.
%
%   This is an internal routine called by DIFFMATRIX_.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('diffmatrixsquare_',xs,ws);
m = length(xs)-1;
D = nan(m+1);
for i = 1:m+1
  for j = [1:i-1 i+1:m+1]
    D(i,j) = ws(j)./(ws(i)*(xs(i)-xs(j)));
  end
  D(i,i) = -sum(D(i,[1:i-1 i+1:end]));
end

