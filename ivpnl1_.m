function [ps,qs,k] = ivpnl1_(f,dfdu,ps,qs,ss,J,K,E1,E2,Ea1,Ea2 ...
                                                                 ,kmax,tol)
%IVPNL1_   Computational kernel for first-order nonlinear IVPs.
%
%   [PS,QS,K] = IVPNL1_(F,DFDU,PS,QS,SS,J,K,E1,E2,EA1,EA2,KMAX,TOL) solves
%
%     U'(T) = F(T,U(T)), U(A) = UA.
%
%   Inputs:
%
%     F: right-hand side of the differential equation
%
%     DFDU: derivative of F(T,U) with respect to U
%
%     PS: sample of a starting guess PSTART(T) on a grid TS of degree M+1
%
%     QS: sample of P' on a grid TS_ of degree M
%
%     SS: degree-M grid of collocation nodes
%
%     J: net-change matrix for the grid TS
%
%     K: truncated integration matrix for the grid TS_
%
%     E1: matrix for resampling from TS to SS
%
%     E2: matrix for resampling from TS_ to SS
%
%     EA1: matrix for evaluating at T = A from a sample on TS
%
%     EA2: matrix for evaluating at T = A from a sample on TS_
%
%     KMAX: maximum number of Newton steps
%
%     TOL: tolerance for termination criterion
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on TS
%
%     QS: sample of the collocation solution derivative P'(T) on TS_
%
%     K: number of Newton steps before convergence, or [] in case of
%     failure
%
%   This is an internal routine called by IVPNL1GEN, IVPNL1UNI, and
%   IVPNL1CHEB.
%
%   Copyright 2019 Brian Sutton

m = length(qs)-1;
% iterate
for k = 1:kmax
  ps_ = E1*ps;
  fs = arrayfun(f,ss,ps_);
  % terminate if solution found
  if all(qs==fs), return; end
  % take a Newton step
  dfdus = arrayfun(dfdu,ss,ps_);
  [dps,dqs] = ivpl1_(-dfdus,ones(m+1,1),fs-qs,[1 0],0 ...
                    ,J,K,E1,E2,Ea1,Ea2);
  ps = ps+dps;
  qs = qs+dqs;
  % signal failure on NaN or infinity
  if any(~isfinite(ps))||any(~isfinite(qs))
    k = []; return;
  end
  % terminate on numerical convergence
  de = tol*max([max(abs(ps)) max(abs(qs)) 1]);
  if max([ max(abs(dps)) max(abs(dqs)) ])<=de, return; end
end
% signal failure to converge
k = [];

