function [ps,qs] = ivpl1cheb(al,be,g,ab,c,d,m)
%IVPL1CHEB   Solve a first-order linear initial-value problem on a
%Chebyshev grid.
%
%   [PS,QS] = IVPL1CHEB(AL,BE,G,[A B],[C0 C1],D,M) computes a collocation
%   solution to the IVP
%
%     AL(T)*U(T) + BE(T)*U'(T) = G(T), C0*U(A) + C1*U'(A) = D
%
%   using a degree-M Chebyshev collocation grid on [A,B].
%
%   Inputs:
%
%     AL: coefficient function of U(T) in the differential equation
%
%     BE: coefficient function of U'(T) in the differential equation
%
%     G: right-hand side of the differential equation
%
%     AB = [A B]: endpoints of the problem domain
%
%     C = [C0 C1]: coefficients in the initial condition
%
%     D: right-hand side of the initial condition
%
%     M: degree of the grid of collocation nodes
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on a Chebyshev grid of
%     degree M+1
%
%     QS: sample of the collocation solution derivative P'(T) on a
%     Chebyshev grid of degree M
%
%   P(T) is a degree-(M+1) polynomial that approximates U(T).
%
%   Copyright 2019 Brian Sutton

narginchk(7,7);
natecheck('ivpl1cheb',al,be,g,ab,c,d,m);
a = ab(1); b = ab(2);
als = samplecheb(al,[a b],m);
bes = samplecheb(be,[a b],m);
gs = samplecheb(g,[a b],m);
[J,K,E1,E2,Ea1,Ea2] = ivp1matcheb(m,b-a);
[ps,qs] = ivpl1_(als,bes,gs,c,d,J,K,E1,E2,Ea1,Ea2);

