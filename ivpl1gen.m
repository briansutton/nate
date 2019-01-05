function [ps,qs] = ivpl1gen(al,be,g,a,c,d,ss,ts,ts_)
%IVPL1GEN   Solve a first-order linear initial-value problem by collocation
%on given grids.
%
%   [PS,QS] = IVPL1GEN(AL,BE,G,A,[C0 C1],D,SS,TS,TS_) computes a
%   collocation solution to the IVP
%
%     AL(T)*U(T) + BE(T)*U'(T) = G(T), C0*U(A) + C1*U'(A) = D
%
%   on [A,B].
%
%   Inputs:
%
%     AL: coefficient function of U(T) in the differential equation
%
%     BE: coefficient function of U'(T) in the differential equation
%
%     G: right-hand side of the differential equation
%
%     A: location of initial condition
%
%     C = [C0 C1]: coefficients in the initial condition
%
%     D: right-hand side of the initial condition
%
%     SS: degree-M grid of collocation nodes
%
%     TS: degree-(M+1) interpolation grid
%
%     TS_: degree-M interpolation grid
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on TS
%
%     QS: sample of the collocation solution derivative P'(T) on TS_
%
%   P(T) is a degree-(M+1) polynomial that approximates U(T).
%
%   No node in TS(2:END) may lie at T = A.
%
%   Copyright 2019 Brian Sutton

narginchk(9,9);
natecheck('ivpl1gen',al,be,g,a,c,d,ss,ts,ts_);
% sample coefficient functions and right-hand side of DE
als = arrayfun(al,ss);
bes = arrayfun(be,ss);
gs = arrayfun(g,ss);
% construct integration, resampling, and evaluation matrices
[J,K] = indefinitegen(ts,ts_);
E1 = evalmatrixgen(ss,ts);
E2 = evalmatrixgen(ss,ts_);
Ea1 = evalmatrixgen(a,ts);
Ea2 = evalmatrixgen(a,ts_);
% solve collocation system
[ps,qs] = ivpl1_(als,bes,gs,c,d,J,K,E1,E2,Ea1,Ea2);

