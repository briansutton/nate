function [ps,qs,rs,k] = ivpnl2_(f,dfdu,dfdv,ps,qs,rs ...
                             ,ss,J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3,kmax,tol)
%IVPNL2_   Computational kernel for second-order nonlinear IVPs.
%
%   [PS,QS,RS,K] = IVPNL2_(F,DFDU,DFDV,PS,QS,RS,SS ...
%                    ,J1,J2,K1,K2,E1,E2,E3,EA1,EA2,EA3,KMAX,TOL)
%
%   solves
%
%     U''(T) = F(T,U(T),U'(T)), U(A) = UA, U'(A) = VA.
%
%   Inputs:
%
%     F: right-hand side of the differential equation
%
%     DFDU: derivative of F(T,U,V) with respect to U
%
%     DFDV: derivative of F(T,U,V) with respect to V
%
%     PS: sample of a starting guess PSTART(T) on a grid TS of degree M+2
%
%     QS: sample of PSTART'(T) on a grid TS_ of degree M+1
%
%     RS: sample of PSTART''(T) on a grid TS__ of degree M
%
%     SS: degree-M grid of collocation nodes
%
%     J1: net-change matrix for the grid TS
%
%     J2: net-change matrix for the grid TS_
%
%     K1: truncated integration matrix for the grid TS_
%
%     K2: truncated integration matrix for the grid TS__
%
%     E1: matrix for resampling from TS to SS
%
%     E2: matrix for resampling from TS_ to SS
%
%     E3: matrix for resampling from TS__ to SS
%
%     EA1: matrix for evaluating at T = A from a sample on TS
%
%     EA2: matrix for evaluating at T = A from a sample on TS_
%
%     EA3: matrix for evaluating at T = A from a sample on TS__
%
%     KMAX: maximum number of Newton steps
%
%     TOL: tolerance for termination criterion
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on TS
%
%     QS: sample of P'(T) on TS_
%
%     RS: sample of P''(T) on TS__
%
%     K: number of Newton steps before convergence, or [] in case of
%     failure
%
%   This is an internal routine called by IVPNL2GEN, IVPNL2UNI, and
%   IVPNL2CHEB.
%
%   Copyright 2019 Brian Sutton

m = length(rs)-1;
for k = 1:kmax
  ps__ = E1*ps;
  qs_ = E2*qs;
  fs = arrayfun(f,ss,ps__,qs_);
  % terminate if solution found
  if all(rs==fs), return; end
  % take a Newton step
  dfdus = arrayfun(dfdu,ss,ps__,qs_);
  dfdvs = arrayfun(dfdv,ss,ps__,qs_);
  [dps,dqs,drs] = ivpl2_(-dfdus,-dfdvs,ones(m+1,1),fs-rs ...
      ,[1 0 0; 0 1 0],[0;0],J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3);
  ps = ps+dps;
  qs = qs+dqs;
  rs = rs+drs;
  % signal failure on NaN or infinity
  if any(~isfinite(ps))||any(~isfinite(qs))||any(~isfinite(rs))
    k = []; return;
  end
  % terminate on numerical convergence
  de = tol*max([ max(abs(ps)) max(abs(qs)) max(abs(rs)) 1]);
  if max([ max(abs(dps)) max(abs(dqs)) max(abs(drs)) ])<=de
    return;
  end
end
k = [];

