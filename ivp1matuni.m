function [J,K,E1,E2,Ea1,Ea2] = ivp1matuni(m,l)
%IVP1MATUNI   Matrices for use in solving first-order IVPs on piecewise-
%uniform grids.
%
%   [J,K,E1,E2,EA1,EA2] = IVP1MATUNI(M,L) constructs integration and
%   evaluation matrices for use in piecewise-uniform collocation on first-
%   order problems.
%
%   Inputs:
%
%     M: degree of a collocation grid
%
%     L: length (B-A)/N of a subinterval in a partition of [A,B] into N
%     parts
%
%   Outputs:
%
%     J: the net-change matrix for a uniform grid of degree M+1
%
%     K: the truncated integration matrix for a uniform grid of degree M
%
%     E1: a matrix for resampling from a uniform grid of degree M+1 to one
%     of degree M
%
%     E2: an identity matrix (because no resampling is necessary)
%
%     EA1: a matrix for evaluating a degree-(M+1) polynomial at its left
%     endpoint, given a sample on a uniform grid
%
%     EA2: a matrix for evaluating a degree-M polynomial at its left
%     endpoint, given a sample on a uniform grid
%
%   This routine is used by IVPL1UNI and IVPNL1UNI.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('ivp1matuni',m,l);
[J,K] = indefiniteuni(m); K = l*K;
E1 = resamplematrixuni(m,m+1);
E2 = eye(m+1);
Ea1 = [ 1 zeros(1,m+1) ];
Ea2 = [ 1 zeros(1,m) ];

