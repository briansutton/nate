function [ps,qs,rs] = ivpl2uni(al,be,ga,g,ab,c,d,m,n)
%IVPL2UNI   Solve a second-order linear initial-value problem on a
%piecewise-uniform grid.
%
%   [PS,QS,RS] = IVPL2UNI(AL,BE,GA,G,[A B],C,D,M,N) computes a collocation
%   solution to the IVP
%
%     AL(T)*U(T) + BE(T)*U'(T) + GA(T)*U''(T) = G(T),
%     C10*U(A) + C11*U'(A) + C12*U''(A) = D1,
%     C20*U(A) + C21*U'(A) + C22*U''(A) = D2
%
%   using an M-by-N piecewise-uniform collocation grid on [A,B].
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
%     M: degree of the grid of collocation nodes on each subinterval
%
%     N: number of subintervals
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on an (M+2)-by-N
%     piecewise-uniform grid
%
%     QS: sample of the collocation solution derivative P'(T) on an (M+1)-
%     by-N piecewise-uniform grid
%
%     RS: sample of the collocation solution second derivative P''(T) on an
%     M-by-N piecewise-uniform grid
%
%   P(T) is a piecewise-polynomial function that approximates U(T). Its
%   pieces are of degree M+2.
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
%     m = 3; n = 1000;
%     [ps,qs,rs] = ivpl2uni(al,be,ga,g,[a b],c,d,m,n);
%     p = interpuni(ps,[a b]);
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(9,9);
natecheck('ivpl2uni',al,be,ga,g,ab,c,d,m,n);
a = ab(1);
b = ab(2);
als = sampleuni(al,[a b],m,n);
bes = sampleuni(be,[a b],m,n);
gas = sampleuni(ga,[a b],m,n);
gs = sampleuni(g,[a b],m,n);
[J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3] = ivp2matuni(m,(b-a)/n);
ps = nan(m+3,n);
qs = nan(m+2,n);
rs = nan(m+1,n);
[ps(:,1),qs(:,1),rs(:,1)] = ivpl2_(        ...
    als(:,1),bes(:,1),gas(:,1),gs(:,1),c,d ...
    ,J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3);
for j = 2:n
  [ps(:,j),qs(:,j),rs(:,j)] = ivpl2_(            ...
      als(:,j),bes(:,j),gas(:,j),gs(:,j)         ...
      ,[1 0 0; 0 1 0],[ps(end,j-1); qs(end,j-1)] ...
      ,J1,J2,K1,K2,E1,E2,E3,Ea1,Ea2,Ea3);
end
ps(1,2:end) = ps(end,1:end-1);
qs(1,2:end) = qs(end,1:end-1);
if m>0, rs(1,2:end) = rs(end,1:end-1); end
