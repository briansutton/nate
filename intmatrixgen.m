function K = intmatrixgen(xs,xs_,a)
%INTMATRIXGEN   Integration matrix for given grids.
%
%   K = INTMATRIXGEN(XS,XS_,A) constructs an integration matrix.
%
%   Inputs:
%
%     XS: degree-(M+1) grid for the antiderivative
%
%     XS_: degree-M grid for the integrand
%
%     A: lower limit of integration
%
%   Outputs:
%
%     K: the integration matrix
%
%   Suppose Q is a degree-M polynomial and P is its antiderivative that
%   satisfies P(A) = 0. If QS is a sample of Q on XS_, then PS = K*QS is a
%   sample of P on XS.
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('intmatrixgen',xs,xs_,a);
ws = baryweights(xs);
K = intmatrix_(xs,ws,xs_,a);

