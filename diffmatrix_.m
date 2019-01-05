function D = diffmatrix_(xs_,xs,ws)
%DIFFMATRIX_   Differentiation matrix for given grids.
%
%   D = DIFFMATRIX_(XS_,XS,WS) computes a differentiation matrix.
%
%   Inputs:
%
%     XS_: grid of degree L on which the derivative Q = P' is sampled
%
%     XS: grid of degree M on which P is sampled
%
%     WS: barycentric weights for the grid XS
%
%   Outputs:
%
%     D: the differentiation matrix
%
%   If PS is a sample of a polynomial P on the grid XS, then QS = D*PS is a
%   sample of Q = P' on the grid XS_.
%
%   This is an internal routine called by DIFFMATRIXGEN and INTMATRIXGEN.
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('diffmatrix_',xs_,xs,ws);
Dxx = diffmatrixsquare_(xs,ws);
E = evalmatrix_(xs_,xs,ws);
D = E*Dxx;

