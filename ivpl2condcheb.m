function kappa = ivpl2condcheb(al,be,ga,ab,c,m)
%IVPL2CONDCHEB   Estimate the absolute condition number of a second-order
%collocation problem on a Chebyshev grid.
%
%   KAPPA = IVPL2CONDCHEB(AL,BE,GA,[A B],[ C10 C11 C12; C20 C21 C22 ],M)
%   estimates the absolute condition number of the Chebyshev collocation
%   problem associated with
%
%     AL(T)*U(T) + BE(T)*U'(T) + GA(T)*U''(T) = G(T),
%     C10*U(A) + C11*U'(A) + C12*U''(A) = D1,
%     C20*U(A) + C21*U'(A) + C22*U''(A) = D2.
%
%   Inputs:
%
%     AL: coefficient function of U(T) in the differential equation
%
%     BE: coefficient function of U'(T) in the differential equation
%
%     GA: coefficient function of U''(T) in the differential equation
%
%     AB = [A B]: endpoints of the problem domain
%
%     C = [ C10 C11 C12; C20 C21 C22 ]: matrix of coefficients in the
%     initial conditions
%
%     M: degree of the grid of collocation nodes
%
%   Outputs:
%
%     KAPPA: an estimate for the absolute condition number of the finite
%     collocation system
%
%   Example: A forced and damped oscillator 2*U + 0.1*U' + U'' = SIN(T).
%
%     al = @(t) 2;
%     be = @(t) 0.1;
%     ga = @(t) 1;
%     g = @(t) sin(t);
%     a = 0; b = 90;
%     c = [ 1 0 0; 0 1 0 ];
%     d = [1; 0];
%     m = 100;
%     [ps,qs,rs] = ivpl2cheb(al,be,ga,g,[a b],c,d,m);
%     p = interpcheb(ps,[a b]);
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(6,6);
natecheck('ivpl2condcheb',al,be,ga,ab,c,m);
s = warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
a = ab(1);
b = ab(2);
als = samplecheb(al,[a b],m);
bes = samplecheb(be,[a b],m);
gas = samplecheb(ga,[a b],m);
[J1,K1] = indefinitecheb(m+1);
[J2,K2] = indefinitecheb(m);
E1 = resamplematrixcheb(m,m+2);
E2 = resamplematrixcheb(m,m+1);
E3 = eye(m+1);
Ea1 = eye(1,m+3);
Ea2 = eye(1,m+2);
Ea3 = eye(1,m+1);
A = diag(als)*E1;
B = diag(bes)*E2;
C = diag(gas)*E3;
L = [ J1             -(b-a)/2*K1 zeros(m+2,m+1) ;
      zeros(m+1,m+3) J2          -(b-a)/2*K2    ;
      A              B           C              ;
      c(1,1)*Ea1     c(1,2)*Ea2  c(1,3)*Ea3     ;
      c(2,1)*Ea1     c(2,2)*Ea2  c(2,3)*Ea3     ];
kappa = condest(L')/norm(L,Inf);
warning(s);

