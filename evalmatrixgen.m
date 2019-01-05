function E = evalmatrixgen(xs_,xs)
%EVALMATRIXGEN   Matrix for evaluating a polynomial at specified points.
%
%   E = EVALMATRIXGEN(XS_,XS) constructs a matrix E for evaluating a
%   polynomial.
%
%   Inputs:
%
%     XS_: points where a polynomial P will be evaluated
%
%     XS: grid of degree M on which P is originally sampled
%
%   Outputs:
%
%     E: the evaluation matrix
%
%   If PS is a sample of a degree-M polynomial on the grid XS, then E*PS is
%   a sample of P on XS_.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('evalmatrixgen',xs_,xs);
ws = baryweights(xs);
E = evalmatrix_(xs_,xs,ws);

