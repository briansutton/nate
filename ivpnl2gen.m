function [ps,qs,rs] = ivpnl2gen(f,dfdu,dfdv,a,ua,va,rstart ...
                                                  ,ss,ts,ts_,ts__,kmax,tol)
%IVPNL2GEN   Solve a second-order nonlinear initial-value problem on given
%grids.
%
%   [PS,QS,RS] = IVPNL2GEN(F,DFDU,DFDV,A,UA,VA,RSTART,SS,TS,TS_,TS__)
%   solves the IVP
%
%     U''(T) = F(T,U(T),U'(T)), U(A) = UA, U'(A) = VA
%
%   using given collocation and interpolation grids.
%
%   Inputs:
%
%     F: right-hand side of the differential equation
%
%     DFDU: derivative of F(T,U,V) with respect to U
%
%     DFDV: derivative of F(T,U,V) with respect to V
%
%     A: location of the initial condition
%
%     UA: initial value U(A)
%
%     VA: initial value U'(A)
%
%     RSTART: starting guess for U''(T)
%
%     SS: degree-M grid of collocation nodes
%
%     TS: degree-(M+2) grid of interpolation nodes
%
%     TS_: degree-(M+1) grid of interpolation nodes
%
%     TS__: degree-M grid of interpolation nodes
%
%     KMAX: maximum number of Newton steps (default: KMAX = 40)
%
%     TOL: tolerance in termination criteria (default: TOL = 1e-12)
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on the grid TS
%
%     QS: sample of the collocation solution derivative P'(T) on the grid
%     TS_
%
%     RS: sample of the collocation solution second derivative P''(T) on
%     the grid TS__
%
%   P(T) is a polynomial approximation to the IVP solution U(T).
%
%   Copyright 2019 Brian Sutton

narginchk(11,13);
if nargin<12, kmax = 40; end
if nargin<13, tol = 1e-12; end
natecheck('ivpnl2gen',f,dfdu,dfdv,a,ua,va,rstart,ss,ts,ts_,ts__,kmax,tol);
rs = arrayfun(rstart,ts__);
qs = antiderivgen(va,rs,ts_,ts__,a);
ps = antiderivgen(ua,qs,ts,ts_,a);
[J1,K1] = indefinitegen(ts,ts_);
[J2,K2] = indefinitegen(ts_,ts__);
E1 = evalmatrixgen(ss,ts);
E2 = evalmatrixgen(ss,ts_);
E3 = evalmatrixgen(ss,ts__);
Ea1 = evalmatrixgen(a,ts);
Ea2 = evalmatrixgen(a,ts_);
Ea3 = evalmatrixgen(a,ts__);
[ps,qs,rs] = ivpnl2_(f,dfdu,dfdv,ps,qs,rs,ss,J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3,kmax,tol);

