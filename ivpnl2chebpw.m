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
%   Copyright 2019 Brian Sutton

narginchk(8,10);
if nargin<9, kmax = 40; end
if nargin<10, tol = 1e-12; end
natecheck('ivpnl2chebpw',f,dfdu,dfdv,ab,ua,va,m,n,kmax,tol);
a = ab(1);
b = ab(2);
[as,bs] = partition_([a b],n);
asbs = [as;bs];
ss = gridchebpw(asbs,m);
ts = gridchebpw(asbs,m+2);
ts_ = gridchebpw(asbs,m+1);
[J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3] = ivp2matcheb(m,(b-a)/n);
ps = nan(m+3,1); qs = nan(m+2,1); rs = nan(m+1,1);
uinit = ua; vinit = va;
warned = false;
%h = waitbar(0,'Computing...');
%h = plot(a,ua,'.');%[[temporary]]
for j = 1:n
  winit = f(as(j),uinit,vinit);
  rsstart = winit*ones(m+1,1);
  qsstart = vinit*ones(m+2,1)+winit*(ts_(:,j)-as(j));
  psstart = uinit*ones(m+3,1)+vinit*(ts(:,j)-as(j))+winit/2*(ts(:,j)-as(j)).^2;
  [ps(:,j),qs(:,j),rs(:,j),k] = ivpnl2_(f,dfdu,dfdv,psstart ...
      ,qsstart,rsstart,ss(:,j),J1,J2,K1,K2,E1               ...
      ,E2,E3,Ea1,Ea2,Ea3,kmax,tol);
  if isempty(k)&&~warned
    warning('NATE:nonlinearIvpFailure' ...
           ,'Nonlinear IVP solver failed to converge');
    warned = true;
  end
  %set(h,'xdata',[get(h,'xdata') as(j)+ts']); set(h,'ydata',[get(h,'ydata') ps(:,j)']); drawnow; %[[temporary]]
  %xlim([a b]);
  %pstep = interpcheb(psstep,[as(j) bs(j)]);
  %qstep = interpcheb(qsstep,[as(j) bs(j)]);
  %rstep = interpcheb(rsstep,[as(j) bs(j)]);
  %ps(as(j)<=ts&ts<=bs(j)) = pstep(ts(as(j)<=ts&ts<=bs(j)));
  %qs(as(j)<=ts_&ts_<=bs(j)) = qstep(ts_(as(j)<=ts_&ts_<=bs(j)));
  %rs(as(j)<=ts__&ts__<=bs(j)) = rstep(ts__(as(j)<=ts__&ts__<=bs(j)));
  %I = as(j)<=ts&ts<=bs(j);
  %ps(I) = interp_(as(j)+xs,ws,pspw,ts(I));
  %I_ = as(j)<=ts_&ts_<=bs(j);
  %qs(I_) = interp_(as(j)+xs_,ws_,qspw,ts_(I_));
  %I__ = as(j)<=ts__&ts__<=bs(j);
  %rs(I__) = interp_(as(j)+xs__,ws__,rspw,ts__(I__));
  uinit = ps(end);
  vinit = qs(end);
  %if mod(j,round(n/10))==0
  %  waitbar(j/n,h);
  %end
end
ps(1,2:end) = ps(end,1:end-1);          % \% continuity @
qs(1,2:end) = qs(end,1:end-1);
if m>0, rs(1,2:end) = rs(end,1:end-1); end

