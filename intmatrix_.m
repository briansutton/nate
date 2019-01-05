function K = intmatrix_(xs,ws,xs_,a)
%INTMATRIX_   Integration matrix on given grids.
%
%   K = INTMATRIX_(XS,WS,XS_,A) constructs an integration matrix on given
%   grids.
%
%   Inputs:
%
%     XS: degree-(M+1) grid on which the antiderivative is sampled
%
%     WS: barycentric weights associated with XS
%
%     XS_: degree-M grid on which the integrand is sampled
%
%     A: lower limit of integration
%
%   Outputs:
%
%     K: the integration matrix
%
%   Suppose Q is a degree-M polynomial and P is its antiderivative that
%   satisfies P(A) = 0. Let QS be a sample of Q on XS_, and let PS be a
%   sample of P on XS. Then PS = K*QS.
%
%   This is an internal routine called by INTMATRIXGEN, INTMATRIXUNI, and
%   INTMATRIXCHEB.
%
%   Copyright 2019 Brian Sutton

narginchk(4,4);
natecheck('intmatrix_',xs,ws,xs_,a);
m = length(xs_)-1;
D = diffmatrix_(xs_,xs,ws);
Ea = evalmatrix_(a,xs,ws);
I = eye(m+1);
K = [ D; Ea ]\[ I; zeros(1,m+1) ];
K(xs==a,:) = 0;

