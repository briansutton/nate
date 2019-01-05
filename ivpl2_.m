function [ps,qs,rs] = ivpl2_(als,bes,gas,gs,c,d ...
                                         ,J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3)
%IVPL2_   Computational kernel for second-order linear IVPs.
%
%   [PS,QS,RS] =
%   IVPL2_(ALS,BES,GAS,GS,C,D,J1,J2,K1,K2,E1,E2,E3,EA1,EA2,EA3) computes a
%   collocation solution to the IVP
%
%     AL(T)*U(T) + BE(T)*U'(T) + GA(T)*U''(T) = G(T),
%     C10*U(A) + C11*U'(A) + C12*U''(A) = D1,
%     C20*U(A) + C21*U'(A) + C22*U''(A) = D2.
%
%   Below, SS denotes a degree-M grid of collocation nodes, and TS, TS_,
%   and TS__ denote interpolation nodes of degrees M+2, M+1, and M,
%   respectively.
%
%   Inputs:
%
%     ALS: sample of the coefficient function AL(T) on SS
%
%     BES: sample of the coefficient function BE(T) on SS
%
%     GAS: sample of the coefficient function GA(T) on SS
%
%     GS: sample of the right-hand side function G(T) on SS
%
%     C = [ C10 C11 C12; C20 C21 C22 ]: matrix of coefficients in the
%     initial conditions
%
%     D = [D1; D2]: vector of right-hand sides for the initial conditions
%
%     J1: net-change matrix for TS
%
%     J2: net-change matrix for TS_
%
%     K1: truncated integration matrix for TS_
%
%     K2: truncated integration matrix for TS__
%
%     E1: matrix for resampling from TS to SS
%
%     E2: matrix for resampling from TS_ to SS
%
%     E3: matrix for resampling from TS__ to SS
%
%     EA1: matrix for evaluating P(A) from a sample PS on TS
%
%     EA2: matrix for evaluating P'(A) from a sample QS on TS_
%
%     EA3: matrix for evaluating P''(A) from a sample RS on TS__
%
%   Outputs
%
%     PS: sample of the collocation solution P(T) on TS
%
%     QS: sample of the collocation solution derivative P'(T) on TS_
%
%     RS: sample of the collocation solution second derivative P''(T) on
%     TS__
%
%   This is an internal routine called by IVPL2GEN, IVPL2UNI, IVPL2CHEB,
%   and IVPNL2_.
%
%   Copyright 2019 Brian Sutton

m = length(als)-1;
A = diag(als)*E1; B = diag(bes)*E2; C = diag(gas)*E3;
L = [ J1             -K1        zeros(m+2,m+1) ;
      zeros(m+1,m+3) J2         -K2            ;
      A              B          C              ;
      c(1,1)*Ea1     c(1,2)*Ea2 c(1,3)*Ea3     ;
      c(2,1)*Ea1     c(2,2)*Ea2 c(2,3)*Ea3     ];
rhs = [ zeros(m+2,1); zeros(m+1,1); gs; d(1); d(2) ];
sol = L\rhs;
ps = sol(1:m+3); qs = sol(m+4:2*m+5); rs = sol(2*m+6:end);

