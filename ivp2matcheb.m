function [J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3] = ivp2matcheb(m,l)
%IVP2MATCHEB   Matrices for use in solving second-order IVPs on Chebyshev
%grids.
%
%   [J1,J2,K1,K2,E1,E2,E3,EA1,EA2,EA3] = IVP2MATCHEB(M,L) constructs
%   integration and evaluation matrices for use in Chebyshev collocation on
%   second-order problems.
%
%   Inputs:
%
%     M: degree of a collocation grid
%
%     L: length B-A of a problem domain [A,B]
%
%   Outputs:
%
%     J1: the net-change matrix for a Chebyshev grid of degree M+2
%
%     K1: the truncated integration matrix for a Chebyshev grid of degree
%     M+1
%
%     J2: the net-change matrix for a Chebyshev grid of degree M+1
%
%     K2: the truncated integration matrix for a Chebyshev grid of degree M
%
%     E1: a matrix for resampling from a Chebyshev grid of degree M+2 to
%     one of degree M
%
%     E2: a matrix for resampling from a Chebyshev grid of degree M+1 to
%     one of degree M
%
%     E3: an identity matrix (because no resampling is necessary)
%
%     EA1: a matrix for evaluating a degree-(M+2) polynomial at its left
%     endpoint, given a sample on a Chebyshev grid
%
%     EA2: a matrix for evaluating a degree-(M+1) polynomial at its left
%     endpoint, given a sample on a Chebyshev grid
%
%     EA3: a matrix for evaluating a degree-M polynomial at its left
%     endpoint, given a sample on a Chebyshev grid
%
%   This routine is used by IVPL2CHEB and IVPNL2CHEBPW.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('ivp2matcheb',m,l);
[J1,K1] = indefinitecheb(m+1);
[J2,K2] = indefinitecheb(m);
K1 = l/2*K1;
K2 = l/2*K2;
E1 = resamplematrixcheb(m,m+2);
E2 = resamplematrixcheb(m,m+1);
E3 = eye(m+1);
Ea1 = [ 1 zeros(1,m+2) ];
Ea2 = [ 1 zeros(1,m+1) ];
Ea3 = [ 1 zeros(1,m) ];

