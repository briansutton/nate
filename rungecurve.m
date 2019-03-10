function [x,y] = rungecurve(ab)
%RUNGECURVE   Plot the Runge curve for a uniform grid.
%
%   RUNGECURVE([A B]) plots the Runge curve in the complex plane for the
%   interval [A,B].
%
%   Inputs:
%
%     AB = [A B]: endpoints of an interval
%
%   If a function has a singularity inside the Runge curve, then the Runge
%   phenomenon is expected.
%
%   Example:
%
%     newfig;
%     rungecurve([-1 1]);
%     axis([-1.5 1.5 -1.5 1.5]);
%     axis square;
%     grid on;
%
%   Copyright 2019 Brian Sutton

narginchk(1,1);
natecheck('rungecurve',ab);
t = linspace(-0.01,0.01,201);
x1 = 1+pi/2*real(-expint(-log(abs(t))));
x1(t==0) = 1;
y1 = t;
f = @(t,u) [2*(atan((1-u(1))/u(2))+atan((1+u(1))/u(2)));
            log((1-u(1))^2+u(2)^2)-log((-1-u(1))^2+u(2)^2)];
sol = ode45(f,[0 0.234],[0;0.5255]);
x2 = deval(sol,linspace(0.234,0,100),1);
y2 = deval(sol,linspace(0.234,0,100),2);
x3 = -x2(end:-1:1);
y3 = y2(end:-1:1);
x4 = -x1;
y4 = -y1;
x5 = -x2;
y5 = -y2;
x6 = -x3;
y6 = -y3;
x7 = x1(1);
y7 = y1(1);
x = [x1 x2 x3 x4 x5 x6 x7];
y = [y1 y2 y3 y4 y5 y6 y7];
a = ab(1);
b = ab(2);
plot((a+b)/2+(b-a)/2*x,(b-a)/2*y,'k');

