function K = intmatrixcheb(m)
%INTMATRIXCHEB   Integration matrix for a Chebyshev grid.
%
%   K = INTMATRIXCHEB(M) constructs an integration matrix for a degree-M
%   Chebyshev grid.
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
%   satisfies P(A) = 0. If QS is a Chebyshev sample of Q on [A,B], then PS
%   = (B-A)/2*K*QS is a Chebyshev sample of P on the same interval.
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('intmatrixcheb',m);
persistent Kstore
if isempty(Kstore), Kstore = cell(1,400); end
if m+1<=length(Kstore)&&~isempty(Kstore{m+1})
  % retrieve matrix if previously computed
  K = Kstore{m+1};
else
  % compute integration matrix
  [xs,ws] = gridcheb([-1 1],m+1);
  xs_ = gridcheb([-1 1],m);
  K = intmatrix_(xs,ws,xs_,-1);
  % store for future use
  if m+1<=length(Kstore)
    Kstore{m+1} = K;
  end
end

