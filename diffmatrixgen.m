function D = diffmatrixgen(xs_,xs)
%DIFFMATRIXGEN   Differentiation matrix for given grids.
%
%   D = DIFFMATRIXGEN(XS_,XS) computes a differentiation matrix.
%
%   Inputs:
%
%     XS_: grid of degree L on which the derivative Q = P' is sampled
%
%     XS: grid of degree M on which P is sampled
%
%   Outputs:
%
%     D: the differentiation matrix
%
%   If PS is a sample of a degree-M polynomial P on the grid XS, then QS =
%   D*PS is a sample of Q = P' on the grid XS_.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('diffmatrixgen',xs_,xs);
ws = baryweights(xs);
D = diffmatrix_(xs_,xs,ws);

