function ws = baryweights(xs)
%BARYWEIGHTS   Barycentric weights for a given grid.
%
%   WS = BARYWEIGHTS(XS) computes a vector of barycentric weights for the
%   grid XS.
%
%   Inputs:
%
%     XS: a vector of distinct nodes on the real line
%
%   Outputs:
%
%     WS: a vector of barycentric weights for use in the barycentric
%     interpolation formula
%
%   Copyright 2020 Brian Sutton

narginchk(1,1);
natecheck('baryweights',xs);
m = length(xs)-1;
C = (max(xs)-min(xs))/4;
ws = nan(m+1,1);
for i = 0:m
  k = [0:i-1 i+1:m];
  ws(i+1) = 1/prod((xs(i+1)-xs(k+1))/C);
end
