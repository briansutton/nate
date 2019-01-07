function bernstein(ab,rho)
%BERNSTEIN   Draw a Bernstein ellipse.
%
%   BERNSTEIN([A B],RHO) draws a Bernstein ellipse on the current axes with
%   foci at Z = A and Z = B and ellipse parameter RHO.
%
%   Inputs:
%
%     AB = [A B]: a vector containing two distinct real numbers for the
%     foci of the ellipse
%
%     RHO: a real number >=1 determining the shape of the ellipse
%
%   RHO = 1 produces a line segment between the foci. Larger values of RHO
%   produce larger and rounder ellipses.
%
%   Example:
%
%     newfig;
%     legend hide;
%     axis([-1.5 1.5 -1.5 1.5]);
%     axis square;
%     bernstein([-1 1],1.2);
%     bernstein([-1 1],1.6);
%     bernstein([-1 1],2.2);
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('bernstein',ab,rho);
a = ab(1);
b = ab(2);
c = (rho^2+1)/(2*rho);
d = rho-c;
t = linspace(0,2*pi,100);
washeld = ishold;
plot((a+b)/2+(b-a)/2*c*cos(t),(b-a)/2*d*sin(t),'k');
hold on;
plot([a b],[0 0],'k');
if ~washeld, hold off; end

