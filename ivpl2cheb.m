function [ps,qs,rs] = ivpl2cheb(al,be,ga,g,ab,c,d,m)
%IVPL2CHEB   Solve a second-order linear initial-value problem on a
%Chebyshev grid.
%
%   [PS,QS,RS] = IVPL2CHEB(AL,BE,GA,G,[A B],C,D,M) computes a collocation
%   solution to
%
%     AL(T)*U(T) + BE(T)*U'(T) + GA(T)*U''(T) = G(T),
%     C10*U(A) + C11*U'(A) + C12*U''(A) = D1,
%     C20*U(A) + C21*U'(A) + C22*U''(A) = D2
%
%   using a degree-M Chebyshev collocation grid on [A,B].
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
%     AB = [A B]: endpoints of the problem domain
%
%     C = [ C10 C11 C12; C20 C21 C22 ]: matrix of coefficients in the
%     initial conditions
%
%     D = [D1; D2]: vector of right-hand sides in the initial conditions
%
%     M: degree of the grid of collocation nodes
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on a Chebyshev grid of
%     degree M+2
%
%     QS: sample of the collocation solution derivative P'(T) on a
%     Chebyshev grid of degree M+1
%
%     RS: sample of the collocation solution second derivative P''(T) on a
%     Chebyshev grid of degree M
%
%   P(T) is a degree-M polynomial that approximates U(T).
%
%   Example: The forced and damped oscillator 2*U + 0.1*U' + U'' = SIN(T).
%
%     al = @(t) 2;
%     be = @(t) 0.1;
%     ga = @(t) 1;
%     g = @(t) sin(t);
%     a = 0; b = 90;
%     c = [ 1 0 0; 0 1 0 ];
%     d = [1; 0];
%     m = 100;
%     [ps,qs,rs] = ivpl2cheb(al,be,ga,g,[a b],c,d,m,n);
%     p = interpcheb(ps,[a b]);
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(8,8);
natecheck('ivpl2cheb',al,be,ga,g,ab,c,d,m);
a = ab(1);
b = ab(2);
als = samplecheb(al,[a b],m);
bes = samplecheb(be,[a b],m);
gas = samplecheb(ga,[a b],m);
gs = samplecheb(g,[a b],m);
[J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3] = ivp2matcheb(m,b-a);
[ps,qs,rs] = ivpl2_(als,bes,gas,gs,c,d,J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3);
