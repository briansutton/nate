function [ps,qs,rs] = ivpl2gen(al,be,ga,g,a,c,d,ss,ts,ts_,ts__)
%IVPL2GEN   Solve a second-order linear initial-value problem by
%collocation on given grids.
%
%   [PS,QS,RS] = IVPL2GEN(AL,BE,GA,G,A,C,D,SS,TS,TS_) computes a
%   collocation solution to the IVP
%
%      AL(T)*U(T) + BE(T)*U'(T) + GA(T)*U''(T) = G(T),
%      C10*U(A) + C11*U'(A) + C12*U''(A) = D1,
%      C20*U(A) + C21*U'(A) + C22*U''(A) = D2.
%
%   on [A,B].
%
%   Inputs:
%
%     AL: coefficient function of U(T) in the differential equation
%
%     BE: coefficient function of U'(T) in the differential equation
%
%     GA: coefficient function of U''(T) in the differential equation
%
%     G: right-hand side of the differential equation
%
%     A: location of initial condition
%
%     C = [ C10 C11 C12; C20 C21 C22 ]: coefficients in the initial
%     conditions
%
%     D = [D1; D2]: right-hand sides for the initial conditions
%
%     SS: degree-M grid of collocation nodes
%
%     TS: degree-(M+2) interpolation grid
%
%     TS_: degree-(M+1) interpolation grid
%
%     TS__: degree-M interpolation grid
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on TS
%
%     QS: sample of the derivative P'(T) on TS_
%
%     RS: sample of the derivative P''(T) on TS__.
%
%   P(T) is a degree-(M+2) polynomial that approximates U(T).
%
%   No node in TS(2:END) or TS_(2:END) may lie at T = A.
%
%   Copyright 2019 Brian Sutton

narginchk(11,11);
natecheck('ivpl2gen',al,be,ga,g,a,c,d,ss,ts,ts_,ts__);
als = arrayfun(al,ss); bes = arrayfun(be,ss); gas = arrayfun(ga,ss);
gs = arrayfun(g,ss);
[J1,K1] = indefinitegen(ts,ts_);
[J2,K2] = indefinitegen(ts_,ts__);
E1 = evalmatrixgen(ss,ts);
E2 = evalmatrixgen(ss,ts_);
E3 = evalmatrixgen(ss,ts__);
Ea1 = evalmatrixgen(a,ts);
Ea2 = evalmatrixgen(a,ts_);
Ea3 = evalmatrixgen(a,ts__);
[ps,qs,rs] = ivpl2_(als,bes,gas,gs,c,d,J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3);

