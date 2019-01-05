function [ps,qs] = ivpnl1gen(f,dfdu,a,ua,qstart,ss,ts,ts_,kmax,tol)
%IVPNL1GEN   Solve a first-order nonlinear initial-value problem on given
%grids.
%
%   [PS,QS] = IVPNL1GEN(F,DFDU,A,UA,QSTART,SS,TS,TS_,KMAX,TOL) solves the
%   IVP
%
%     U'(T) = F(T,U(T)), U(A) = UA
%
%   using given collocation and interpolation grids.
%
%   Inputs:
%
%     F: right-hand side of the differential equation
%
%     DFDU: derivative of F(T,U) with respect to U
%
%     A: location of the initial condition
%
%     UA: initial value
%
%     QSTART: starting guess for U'(T)
%
%     SS: degree-M grid of collocation nodes
%
%     TS: degree-(M+1) grid of interpolation nodes
%
%     TS_: degree-M grid of interpolation nodes
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
%   P(T) is a polynomial approximation to the IVP solution U(T).
%
%   Copyright 2019 Brian Sutton

narginchk(8,10);
if nargin<9, kmax = 40; end
if nargin<10, tol = 1e-12; end
natecheck('ivpnl1gen',f,dfdu,a,ua,qstart,ss,ts,ts_,kmax,tol);
% form starting guess
qs = arrayfun(qstart,ts_);
ps = antiderivgen(ua,qs,ts,ts_,a);
% construct integration, resampling, and evaluation matrices
[J,K] = indefinitegen(ts,ts_);
E1 = evalmatrixgen(ss,ts);
E2 = evalmatrixgen(ss,ts_);
Ea1 = evalmatrixgen(a,ts);
Ea2 = evalmatrixgen(a,ts_);
% solve by collocation and Newton's method
[ps,qs,k] = ivpnl1_(f,dfdu,ps,qs,ss,J,K,E1,E2,Ea1,Ea2,kmax,tol);
if isempty(k)
  warning('NATE:nonlinearIvpFailure' ...
         ,'Nonlinear IVP solver failed to converge');
end

