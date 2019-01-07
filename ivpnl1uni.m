function [ps,qs] = ivpnl1uni(f,dfdu,ab,ua,m,n,kmax,tol)
%IVPNL1UNI   Solve a first-order nonlinear IVP on a piecewise-uniform grid.
%
%   [PS,QS] = IVPNL1UNI(F,DFDU,[A B],UA,M,N,KMAX,TOL) solves the IVP
%
%     U'(T) = F(T,U(T)), U(A) = UA
%
%   using piecewise-uniform grids.
%
%   Inputs:
%
%     F: right-hand side of the differential equation
%
%     DFDU: derivative of F(T,U) with respect to U
%
%     AB = [A B]: endpoints of the problem domain
%
%     UA: initial value
%
%     M: degree of the grid of collocation nodes on each subinterval
%
%     N: number of subintervals
%
%     KMAX: maximum number of Newton steps on each subinterval (default:
%     KMAX = 10)
%
%     TOL: tolerance in termination criteria (default: TOL = 1e-12)
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on a piecewise-uniform
%     grid
%
%     QS: sample of the collocation solution derivative P'(T) on a
%     piecewise-uniform grid
%
%   P(T) is a piecewise-polynomial approximation to the IVP solution U(T).
%   It is sampled on a degree-(M+1) grid on each subinterval to produce PS.
%   Its derivative is sampled on a degree-M grid on each subinterval to
%   produce QS.
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
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(6,8);
if nargin<7, kmax = 10; end
if nargin<8, tol = 1e-12; end
natecheck('ivpnl1uni',f,dfdu,ab,ua,m,n,kmax,tol);
a = ab(1); b = ab(2); [as,~] = partition_([a b],n);
ss = griduni([a b],m,n);
xs = griduni([0 (b-a)/n],m+1);
[J,K,E1,E2,Ea1,Ea2] = ivp1matuni(m,(b-a)/n);
ps = nan(m+2,n); qs = nan(m+1,n);
warned = false;
% set initial condition for first subinterval
uinit = ua;
% march across subintervals
for j = 1:n
  % construct starting guess on subinterval
  vinit = f(as(j),uinit);
  qsstart = repmat(vinit,m+1,1);
  psstart = repmat(uinit,m+2,1)+vinit*xs;
  % solve on subinterval
  [ps(:,j),qs(:,j),k] = ivpnl1_(f,dfdu,psstart,qsstart,ss(:,j) ...
                               ,J,K,E1,E2,Ea1,Ea2,kmax,tol);
  if isempty(k)&&~warned
    warning('NATE:nonlinearIvpFailure' ...
           ,'Nonlinear IVP solver failed to converge');
    warned = true;
  end
  % set initial condition on next subinterval
  uinit = ps(end,j);
end
% enforce continuity
ps(1,2:end) = ps(end,1:end-1);
if m>0, qs(1,2:end) = qs(end,1:end-1); end

