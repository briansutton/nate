function [ps,qs,asbs] = ivpnl1chebpw(f,dfdu,ab,ua,m,n,kmax,tol)
%IVPNL1CHEBPW   Solve a first-order nonlinear IVP on a piecewise-Chebyshev
%grid.
%
%   [PS,QS,ASBS] = IVPNL1CHEBPW(F,DFDU,[A B],UA,M,N,KMAX,TOL) solves the
%   IVP
%
%     U'(T) = F(T,U(T)), U(A) = UA
%
%   using piecewise-Chebyshev grids.
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
%     KMAX = 40)
%
%     TOL: tolerance in termination criteria (default: TOL = 1e-12)
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on a piecewise-Chebyshev
%     grid
%
%     QS: sample of the collocation solution derivative P'(T) on a
%     piecewise-Chebyshev grid
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     in the piecewise-defined grid
%
%   The subintervals are all of the same width.
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
%     m = 32; n = 10;
%     [ps,qs,asbs] = ivpnl1chebpw(f,dfdu,[a b],ua,m,n);
%     p = interpchebpw(ps,asbs);
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(6,8);
if nargin<7, kmax = 40; end
if nargin<8, tol = 1e-12; end
natecheck('ivpnl1chebpw',f,dfdu,ab,ua,m,n,kmax,tol);
a = ab(1); b = ab(2);
[as,bs] = partition_([a b],n); asbs = [as;bs];
ss = gridchebpw(asbs,m);
ts = gridchebpw(asbs,m+1);
[J,K,E1,E2,Ea1,Ea2] = ivp1matcheb(m,(b-a)/n);
ps = nan(m+2,n); qs = nan(m+1,n);
warned = false;
% set initial condition for first subinterval
uinit = ua;
% march across subintervals
for j = 1:n
  % construct starting guess on subinterval
  vinit = f(as(j),uinit);
  qsstart = repmat(vinit,m+1,1);
  psstart = repmat(uinit,m+2,1)+vinit*(ts(:,j)-as(j));
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
qs(1,2:end) = qs(end,1:end-1);

