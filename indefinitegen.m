function [J,K] = indefinitegen(xs,xs_)
%INDEFINITEGEN   Indefinite integral on given grids.
%
%   [J,K] = INDEFINITEGEN(XS,XS_) computes a pair of matrices J, K for
%   testing if one polynomial is an antiderivative of another.
%
%   Inputs:
%
%     XS: a grid of degree M+1
%
%     XS_: a grid of degree M
%
%   Outputs:
%
%     J: matrix for evaluating P(X) - P(A) on the grid XS
%
%     K: matrix for evaluating the integral of Q(X) on the grid XS
%
%   If PS is a sample of a degree-(M+1) polynomial P on XS and QS is a
%   sample of a degree-M polynomial Q on XS_, then P is an antiderivative
%   of Q if and only if J*PS = K*QS.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('indefinitegen',xs,xs_);
m = length(xs_)-1;
J = [ -ones(m+1,1) eye(m+1) ];
K = intmatrixgen(xs,xs_,xs(1));
K = K(2:end,:);

