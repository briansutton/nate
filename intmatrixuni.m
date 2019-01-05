function K = intmatrixuni(m)
%INTMATRIXUNI   Integration matrix for a uniform grid.
%
%   K = INTMATRIXUNI(M) constructs an integration matrix for a degree-M
%   uniform grid.
%
%   Inputs:
%
%     M: degree of a grid
%
%   Outputs:
%
%     K: the integration matrix
%
%   Suppose Q is a degree-M polynomial and P is its antiderivative that
%   satisfies P(A) = 0. If QS is a sample of Q on a uniform grid over
%   [A,B], then PS = (B-A)*K*QS is a sample of P on the same interval.
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('intmatrixuni',m);
persistent Kstore
if isempty(Kstore), Kstore = cell(1,9); end
if m+1<=length(Kstore)&&~isempty(Kstore{m+1})
  % retrieve matrix if previously computed
  K = Kstore{m+1};
else
  % compute integration matrix
  [xs,ws] = griduni([0 1],m+1,1);
  xs_ = griduni([0 1],m,1);
  K = intmatrix_(xs,ws,xs_,0);
  % store for future use
  if m+1<=length(Kstore)
    Kstore{m+1} = K;
  end
end

