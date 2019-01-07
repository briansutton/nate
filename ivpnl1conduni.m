function kappa = ivpnl1conduni(dfdu,ab,m,n,p)
%IVPNL1CONDUNI   Absolute condition number for a first-order nonlinear IVP
%on a piecewise-uniform grid.
%
%   KAPPA = IVPNL1CONDUNI(DFDU,[A B],M,N,P) estimates the absolute
%   condition number for a piecewise-uniform collocation problem that
%   discretizes
%
%     U'(T) = F(T,U(T)), U(A) = UA.
%
%   Inputs:
%
%     DFDU: derivative of the right-hand side F(T,U) of the differential
%     equation with respect to U
%
%     AB = [A B]: endpoints of the problem domain
%
%     M: degree of the grid of collocation nodes on each subinterval
%
%     N: number of subintervals
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
%     m = 3; n = 1000;
%     [ps,qs] = ivpnl1uni(f,dfdu,[a b],ua,m,n);
%     p = interpuni(ps,[a b]);
%     ivpnl1conduni(dfdu,[a b],m,n,p)
%
%   Copyright 2019 Brian Sutton

narginchk(5,5);
natecheck('ivpnl1conduni',dfdu,ab,m,n,p);
al = @(t) -dfdu(t,p(t));
be = @(t) 1;
kappa = ivpl1conduni(al,be,ab,[1 0],m,n);

