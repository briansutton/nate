function [J,K,E1,E2,Ea1,Ea2] = ivp1matcheb(m,l)
%IVP1MATCHEB   Matrices for use in solving first-order IVPs on Chebyshev
%grids.
%
%   [J,K,E1,E2,EA1,EA2] = IVP1MATCHEB(M,L) constructs integration and
%   evaluation matrices for use in Chebyshev collocation on first-order
%   problems.
%
%   Inputs:
%
%     M: degree of a collocation grid
%
%     L: length B-A of a problem domain [A,B]
%
%   Outputs:
%
%     J: the net-change matrix for a Chebyshev grid of degree M+1
%
%     K: the truncated integration matrix for a Chebyshev grid of degree M
%
%     E1: a matrix for resampling from a Chebyshev grid of degree M+1 to
%     one of degree M
%
%     E2: an identity matrix (because no resampling is necessary)
%
%     EA1: a matrix for evaluating a degree-(M+1) polynomial at its left
%     endpoint, given a sample on a Chebyshev grid
%
%     EA2: a matrix for evaluating a degree-M polynomial at its left
%     endpoint, given a sample on a Chebyshev grid
%
%   This routine is used by IVPL1CHEB and IVPNL1CHEBPW.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('ivp1matcheb',m,l);
[J,K] = indefinitecheb(m); K = l/2*K;
E1 = resamplematrixcheb(m,m+1);
E2 = eye(m+1);
Ea1 = [ 1 zeros(1,m+1) ];
Ea2 = [ 1 zeros(1,m) ];

