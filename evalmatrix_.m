function E = evalmatrix_(xs_,xs,ws)
%EVALMATRIX_   Evaluation matrix.
%
%   E = EVALMATRIX_(XS_,XS,WS) constructs a matrix E for evaluating a
%   polynomial.
%
%   Inputs:
%
%     XS_: points where a polynomial P will be evaluated
%
%     XS: grid of degree M on which P is originally sampled
%
%     WS: barycentric weights for the grid XS
%
%   Outputs:
%
%     E: the evaluation matrix
%
%   If PS is a sample of a degree-M polynomial on the grid XS, then E*PS is
%   a sample of P on XS_.
%
%   This is an internal routine called by EVALMATRIXGEN, RESAMPLEMATRIXUNI,
%   and RESAMPLEMATRIXCHEB.
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('evalmatrix_',xs_,xs,ws);
m = length(xs_)-1;
l = length(xs)-1;
E = nan(m+1,l+1);
for i = 1:m+1
  if any(xs==xs_(i))
    E(i,:) = xs_(i)==xs;
  else
    E(i,:) = ws./(xs_(i)-xs);
    E(i,:) = E(i,:)/sum(E(i,:));
  end
end

