function kappa = ivpl1condcheb(al,be,ab,c,m)
%IVPL1CONDCHEB   Estimate the absolute condition number of a first-order
%collocation problem on a Chebyshev grid.
%
%   KAPPA = IVPL1CONDCHEB(AL,BE,[A B],[C0 C1],M) estimates the absolute
%   condition number of the Chebyshev collocation problem associated with
%
%     AL(T)*U(T) + BE(T)*U'(T) = G(T), C0*U(A) + C1*U'(A) = D.
%
%   Inputs:
%
%     AL: coefficient function of U(T) in the differential equation
%
%     BE: coefficient function of U'(T) in the differential equation
%
%     AB = [A B]: endpoints of the problem domain
%
%     C = [C0 C1]: vector of coefficients in the initial condition
%
%     M: degree of the grid of collocation nodes
%
%   Outputs:
%
%     KAPPA: an estimate for the absolute condition number of the finite
%     collocation system
%
%   Example: T*U + U' = G(T).
%
%     al = @(t) t;
%     be = @(t) 1;
%     a = 0; b = 6;
%     c = [1 0];
%     m = 30;
%     ivpl1condcheb(al,be,[a b],c,m)
%
%   Copyright 2019 Brian Sutton

narginchk(5,5);
natecheck('ivpl1condcheb',al,be,ab,c,m);
s = warning('query','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
a = ab(1); b = ab(2);
als = samplecheb(al,[a b],m);
bes = samplecheb(be,[a b],m);
[J,K,E1,E2,Ea1,Ea2] = ivp1matcheb(m,b-a);
A = diag(als)*E1; B = diag(bes)*E2;
c0 = c(1); c1 = c(2);
L = [ J      -K     ;
      A      B      ;
      c0*Ea1 c1*Ea2 ];
kappa = condest(L')/norm(L,Inf);
warning(s);

