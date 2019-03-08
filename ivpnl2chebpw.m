function [ps,qs,rs,asbs] = ivpnl2chebpw(f,dfdu,dfdv,ab,ua,va,m,n,kmax,tol)
%IVPNL2CHEBPW   Solve a second-order nonlinear IVP on a piecewise-Chebyshev
%grid.
%
%   [PS,QS,RS,ASBS] = IVPNL2CHEBPW(F,DFDU,DFDV,[A B],UA,VA,M,N,KMAX,TOL)
%   solves the IVP
%
%     U''(T) = F(T,U(T),U'(T)), U(A) = UA, U'(A) = VA
%
%   using piecewise-Chebyshev grids.
%
%   Inputs:
%
%     F: right-hand side of the differential equation
%
%     DFDU: derivative of F(T,U,V) with respect to U
%
%     DFDV: derivative of F(T,U,V) with respect to V
%
%     AB = [A B]: endpoints of the problem domain
%
%     UA: initial value U(A)
%
%     VA: initial value U'(A)
%
%     M: degree of the grid of collocation nodes on each subinterval
%
%     N: number of subintervals
%
%     KMAX: maximum number of Newton steps (default: KMAX = 40)
%
%     TOL: tolerance in termination criteria (default: TOL = 1e-12)
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on a piecewise-Chebyshev
%     grid
%
%     QS: sample of P'(T) on a piecewise-Chebyshev grid
%
%     QS: sample of P''(T) on a piecewise-Chebyshev grid
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     in the piecewise-defined grid
%
%   The subintervals are all of the same width.
%
%   P(T) is a piecewise-polynomial approximation to the IVP solution U(T).
%   It is sampled on a degree-(M+2) grid on each subinterval to produce PS.
%   Its first and second derivatives are sampled on grids of degrees M+1
%   and M, respectively, to produce QS and RS.
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
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(8,10);
if nargin<9, kmax = 40; end
if nargin<10, tol = 1e-12; end
natecheck('ivpnl2chebpw',f,dfdu,dfdv,ab,ua,va,m,n,kmax,tol);
a = ab(1); b = ab(2);
[as,bs] = partition_([a b],n); asbs = [as;bs];
ss = gridchebpw(asbs,m);
ts = gridchebpw(asbs,m+2);
ts_ = gridchebpw(asbs,m+1);
[J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3] = ivp2matcheb(m,(b-a)/n);
ps = nan(m+3,1); qs = nan(m+2,1); rs = nan(m+1,1);
warned = false;
uinit = ua; vinit = va;
for j = 1:n
  winit = f(as(j),uinit,vinit);
  rsstart = repmat(winit,m+1,1);
  qsstart = repmat(vinit,m+2,1)+winit*(ts_(:,j)-as(j));
  psstart = repmat(uinit,m+3,1)+vinit*(ts(:,j)-as(j))+winit/2*(ts(:,j)-as(j)).^2;
  [ps(:,j),qs(:,j),rs(:,j),k] = ivpnl2_(f,dfdu,dfdv ...
      ,psstart,qsstart,rsstart,ss(:,j),J1,J2,K1,K2  ...
      ,E1,E2,E3,Ea1,Ea2,Ea3,kmax,tol);
  if isempty(k)&&~warned
    warning('NATE:nonlinearIvpFailure' ...
           ,'Nonlinear IVP solver failed to converge');
    warned = true;
  end
  uinit = ps(end,j); vinit = qs(end,j);
end
ps(1,2:end) = ps(end,1:end-1);
qs(1,2:end) = qs(end,1:end-1);
rs(1,2:end) = rs(end,1:end-1);

