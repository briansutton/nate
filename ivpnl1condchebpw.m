function kappa = ivpnl1condchebpw(dfdu,asbs,m,p)
%IVPNL1CONDCHEBPW   Absolute condition number for a first-order nonlinear
%IVP on a piecewise-Chebyshev grid.
%
%   KAPPA = IVPNL1CONDCHEBPW(DFDU,ASBS,M,P) estimates the absolute
%   condition number for a finite collocation system that discretizes
%
%     U'(T) = F(T,U(T)), U(A) = UA.
%
%   Inputs:
%
%     DFDU: derivative of the right-hand side F(T,U) of the differential
%     equation with respect to U
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     in the partition of the domain [A,B]
%
%     M: degree of the grid of collocation nodes on each subinterval
%
%     P: an estimate for the solution U
%
%   Outputs:
%
%     KAPPA: estimated absolute condition number
%
%   Example: The logistic equation U' = U*(1-U).
%
%     f = @(t,u) u*(1-u);
%     dfdu = @(t,u) 1-2*u;
%     a = 0; b = 10;
%     ua = 0.1;
%     m = 32; n = 10;
%     [ps,qs,asbs] = ivpnl1chebpw(f,dfdu,[a b],ua,m,n);
%     p = interpchebpw(ps,asbs);
%     ivpnl1condchebpw(dfdu,asbs,m,p)
%
%   Copyright 2019 Brian Sutton

narginchk(4,4);
natecheck('ivpnl1condchebpw',dfdu,asbs,m,p);
al = @(t) -dfdu(t,p(t));
be = @(t) 1;
kappa = ivpl1condchebpw(al,be,asbs,[1 0],m);

