function kappa = ivpnl2conduni(dfdu,dfdv,ab,m,n,p,q)
%IVPNL2CONDUNI   Absolute condition number for a second-order nonlinear IVP
%on a piecewise-uniform grid.
%
%   KAPPA = IVPNL2CONDUNI(DFDU,DFDV,[A B],M,N,P,Q) estimates the absolute
%   condition number for a finite collocation system that discretizes
%
%     U''(T) = F(T,U(T),U'(T)), U(A) = UA, U'(A) = VA
%
%   on [A,B].
%
%   Inputs:
%
%     DFDU: derivative of the right-hand side F(T,U,V) of the differential
%     equation with respect to U
%
%     DFDV: derivative of the right-hand side F(T,U,V) of the differential
%     equation with respect to V
%
%     AB = [A B]: endpoints of the problem domain [A,B]
%
%     M: degree of the grid of collocation nodes on each subinterval
%
%     N: number of subintervals
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
%     m = 3; n = 1000;
%     [ps,qs,rs] = ivpnl2uni(f,dfdu,dfdv,[a b],ua,va,m,n);
%     p = interpuni(ps,[a b]);
%     q = interpuni(qs,[a b]);
%     ivpnl2conduni(dfdu,dfdv,[a b],m,n,p,q)
%
%   Copyright 2019 Brian Sutton

narginchk(7,7);
natecheck('ivpnl2conduni',dfdu,dfdv,ab,m,n,p,q);
al = @(t) -dfdu(t,p(t),q(t));
be = @(t) -dfdv(t,p(t),q(t));
ga = @(t) 1;
kappa = ivpl2conduni(al,be,ga,ab,[ 1 0 0; 0 1 0 ],m,n);

