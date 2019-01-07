function [ps,qs] = ivpl1uni(al,be,g,ab,c,d,m,n)
%IVPL1UNI   Solve a first-order linear initial-value problem on a
%piecewise-uniform grid.
%
%   [PS,QS] = IVPL1UNI(AL,BE,G,[A B],[C0 C1],D,M,N) computes a collocation
%   solution to the IVP
%
%     AL(T)*U(T) + BE(T)*U'(T) = G(T), C0*U(A) + C1*U'(A) = D
%
%   on [A,B] using piecewise-uniform grids.
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
%     M: degree of the grid of collocation nodes on each subinterval
%
%     N: number of subintervals
%
%   Outputs:
%
%     PS: sample of the collocation solution P(T) on an (M+1)-by-N
%     piecewise-uniform grid
%
%     QS: sample of the collocation solution derivative P'(T) on an M-by-N
%     piecewise-uniform grid
%
%   P(T) is a piecewise-polynomial function that approximates U(T). Its
%   pieces are of degree M+1.
%
%   Example: T*U + U' = 0.
%
%     al = @(t) t;
%     be = @(t) 1;
%     g = @(t) 0;
%     a = 0; b = 6;
%     c = [1 0];
%     d = 1;
%     m = 2; n = 200;
%     [ps,qs] = ivpl1uni(al,be,g,[a b],c,d,m,n);
%     p = interpuni(ps,[a b]);
%     newfig;
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(8,8);
natecheck('ivpl1uni',al,be,g,ab,c,d,m,n);
a = ab(1); b = ab(2);
% sample coefficient functions and right-hand side of DE
als = sampleuni(al,[a b],m,n);
bes = sampleuni(be,[a b],m,n);
gs = sampleuni(g,[a b],m,n);
% construct integration, resampling, and evaluation matrices
[J,K,E1,E2,Ea1,Ea2] = ivp1matuni(m,(b-a)/n);
% solve one subinterval at a time
ps = nan(m+2,n); qs = nan(m+1,n);
[ps(:,1),qs(:,1)] = ivpl1_(als(:,1),bes(:,1) ...
    ,gs(:,1),c,d,J,K,E1,E2,Ea1,Ea2);
for j = 2:n
  [ps(:,j),qs(:,j)] = ivpl1_(als(:,j),bes(:,j),gs(:,j),[1 0] ...
      ,ps(end,j-1),J,K,E1,E2,Ea1,Ea2);
end
% enforce continuity
ps(1,2:end) = ps(end,1:end-1);
if m>0, qs(1,2:end) = qs(end,1:end-1); end

