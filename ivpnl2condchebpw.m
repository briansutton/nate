function kappa = ivpnl2condchebpw(dfdu,dfdv,asbs,m,p,q)
%IVPNL2CONDCHEBPW   Absolute condition number for a second-order nonlinear
%IVP on a piecewise-Chebyshev grid.
%
%   KAPPA = IVPNL2CONDCHEBPW(DFDU,DFDV,ASBS,M,P,Q) estimates the absolute
%   condition number for a finite collocation system that discretizes
%
%     U''(T) = F(T,U(T),U'(T)), U(A) = UA, U'(A) = VA.
%
%   Inputs:
%
%     DFDU: derivative of the right-hand side F(T,U,V) of the differential
%     equation with respect to U
%
%     DFDV: derivative of the right-hand side F(T,U,V) of the differential
%     equation with respect to V
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     in the partition of the domain [A,B]
%
%     M: degree of the grid of collocation nodes on each subinterval
%
%     P: an estimate for the solution U(T)
%
%     Q: an estimate for the solution derivative U'(T)
%
%   Outputs:
%
%     KAPPA: estimated absolute condition number
%
%   Example: The Van der Pol equation U'' = -U + 5*(1-U^2)*U'.
%
%     f = @(t,u,v) -u+5*(1-u^2)*v;
%     dfdu = @(t,u,v) -1-10*u*v;
%     dfdv = @(t,u,v) 5*(1-u^2);
%     a = 0; b = 50;
%     ua = 0.1; va = 0;
%     m = 32; n = 100;
%     [ps,qs,rs,asbs] = ivpnl2chebpw(f,dfdu,dfdv,[a b],ua,va,m,n);
%     p = interpchebpw(ps,asbs);
%     q = interpchebpw(qs,asbs);
%     ivpnl2condchebpw(dfdu,dfdv,asbs,m,p,q)
%
%   Copyright 2019 Brian Sutton

narginchk(6,6);
natecheck('ivpnl2condchebpw',dfdu,dfdv,asbs,m,p,q);
al = @(t) -dfdu(t,p(t),q(t));
be = @(t) -dfdv(t,p(t),q(t));
ga = @(t) 1;
kappa = ivpl2condchebpw(al,be,ga,asbs,[1 0 0; 0 1 0],m);

