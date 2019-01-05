function [ps,qs] = ivpl1_(als,bes,gs,c,d,J,K,E1,E2,Ea1,Ea2)
%IVPL1_   Solve a first-order linear initial-value problem.
%
%   [PS,QS] = IVPL1_(ALS,BES,GS,[C0 C1],D,J,K,E1,E2,EA1,EA2) computes a
%   collocation solution to the IVP
%
%     AL(T)*U(T) + BE(T)*U'(T) = G(T), C0*U(A) + C1*U'(A) = D.
%
%   Below, SS denotes a degree-M grid of collocation nodes, and TS and TS_
%   denote interpolation grids of degrees M+1 and M, respectively.
%
%   Inputs:
%
%     ALS: sample of the coefficient function AL(T) on SS
%
%     BES: sample of the coefficient function BE(T) on SS
%
%     GS: sample of the right-hand side function G(T) on SS
%
%     C = [C0 C1]: coefficients in the initial condition
%
%     D: right-hand side of the initial condition
%
%     J: net-change matrix for TS
%
%     K: truncated integration matrix for TS_
%
%     E1: matrix for resampling from TS to SS
%
%     E2: matrix for resampling from TS_ to SS
%
%     EA1: matrix for evaluating P(A) from a sample PS on TS
%
%     EA2: matrix for evaluating P'(A) from a sample QS on TS_
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on TS
%
%     QS: sample of the collocation solution derivative P'(T) on TS_
%
%   This is an internal routine called by IVPL1GEN, IVPL1UNI, IVLP1CHEB,
%   and IVPNL1_.
%
%   Copyright 2019 Brian Sutton

m = length(als)-1;
A = diag(als)*E1; B = diag(bes)*E2;
c0 = c(1); c1 = c(2);
L = [ J      -K     ;
      A      B      ;
      c0*Ea1 c1*Ea2 ];
rhs = [ zeros(m+1,1); gs; d ];
sol = L\rhs;
ps = sol(1:m+2); qs = sol(m+3:end);

