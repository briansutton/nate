function y = interpuni_(xs,ws,ps,a,b,l,n,x)
%INTERPUNI_   Evaluate an interpolating piecewise-polynomial function.
%
%   Y = INTERPUNI_(XS,WS,PS,A,B,L,N,X) constructs a piecewise-polynomial
%   function P that interpolates points on a piecewise-uniform grid and
%   evaluates P(X).
%
%   Inputs:
%
%     XS: a uniform grid on [0,1]
%
%     WS: barycentric weights for XS
%
%     PS: vertical coordinates of interpolation points
%
%     A: left endpoint of interpolation interval
%
%     B: right endpoint of interpolation interval
%
%     L: width of each subinterval in the piecewise-uniform grid
%
%     N: the number of subintervals in the piecewise-uniform grid
%
%     X: the location where the interpolating P should be evaluated
%
%   Outputs:
%
%     Y: the value Y = P(X) of the interpolating function
%
%   This is an internal routine called by INTERPUNI.
%
%   Copyright 2019 Brian Sutton

y = nan(size(x));
for k = 1:length(x)
  if x(k)>=a&&x(k)<=b
    % find subinterval
    j = min(floor((x(k)-a)/l)+1,n);
    % evaluate polynomial piece
    t = (x(k)-(a+(j-1)*l))/l;
    y(k) = interp_(xs,ws,ps(:,j),t);
  end
end

