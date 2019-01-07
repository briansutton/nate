function [ps,qs,rs] = ivpnl2uni(f,dfdu,dfdv,ab,ua,va,m,n,kmax,tol)
%IVPNL2UNI   Solve a second-order nonlinear initial-value problem on a
%piecewise-uniform grid.
%
%   IVPNL2UNI solves the second-order nonlinear IVP
%
%     U''(T) = F(T,U(T),U'(T)), U(A) = UA, U'(A) = VA
%
%   on an interval [A,B] using piecewise-uniform grids.
%
%   Inputs:
%
%     F: right-hand side of the differential equation
%
%     DFDU: partial derivative of F(T,U,V) with respect to U
%
%     DFDV: partial derivative of F(T,U,V) with respect to V
%
%     AB = [A B]: endpoints of problem domain
%
%     UA: initial value U(A)
%
%     VA: initial value U'(A)
%
%     M: polynomial degree
%
%     N: number of subintervals
%
%     KMAX: maximum number of Newton steps on each subinterval (default:
%     KMAX = 10)
%
%     TOL: tolerance for Newton iteration (default: TOL = 1e-12)
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on an (M+2)-by-N
%     piecewise-uniform grid
%
%     QS: sample of the collocation solution derivative P'(T) on an (M+1)-
%     by-N piecewise-uniform grid
%
%     RS: sample of the collocation solution second derivative P''(T) on an
%     M-by-N piecewise-uniform grid
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
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(8,10);
if nargin<9, kmax = 10; end
if nargin<10, tol = 1e-12; end
natecheck('ivpnl2uni',f,dfdu,dfdv,ab,ua,va,m,n,kmax,tol);
a = ab(1); b = ab(2); [as,~] = partition_([a b],n);
ss = griduni([a b],m,n);
xs = griduni([0 (b-a)/n],m+2); xs_ = griduni([0 (b-a)/n],m+1);
[J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3] = ivp2matuni(m,(b-a)/n);
ps = nan(m+3,n); qs = nan(m+2,n); rs = nan(m+1,n);
uinit = ua; vinit = va;
warned = false;
for j = 1:n
  winit = f(as(j),uinit,vinit);
  rsstart = repmat(winit,m+1,1);
  qsstart = repmat(vinit,m+2,1)+winit*xs_;
  psstart = repmat(uinit,m+3,1)+vinit*xs+winit/2*xs.^2;
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
ps(1,2:end) = ps(end,1:end-1); qs(1,2:end) = qs(end,1:end-1);
if m>0, rs(1,2:end) = rs(end,1:end-1); end

