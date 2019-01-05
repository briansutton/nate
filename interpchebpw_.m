function y = interpchebpw_(xs,ws,ps,as,bs,x)
%INTERPCHEBPW_   Evaluate a piecewise-polynomial interpolant that is
%defined on a piecewise-Chebyshev grid.
%
%   Y = INTERPCHEBPW_(XS,WS,PS,AS,BS,X) constructs a piecewise-polynomial
%   function P that interpolates points on a piecewise-Chebyshev grid and
%   evaluates P(X).
%
%   Inputs:
%
%     XS: a Chebyshev grid on [0,1]
%
%     WS: barycentric weights for XS
%
%     PS: vertical coordinates of interpolation points
%
%     AS: left endpoints of subintervals
%
%     BS: right endpoints of subintervals
%
%     X: the location where the interpolating P should be evaluated
%
%   Outputs:
%
%     Y: the value Y = P(X) of the interpolating function
%
%   This is an internal routine called by INTERPCHEBPW.
%
%   Copyright 2019 Brian Sutton

y = nan(size(x));
for k = 1:length(x)
  if x(k)>=as(1)&&x(k)<=bs(end)
    % find the subinterval
    j = find(as>x(k),1);
    if isempty(j)
      j = length(as);
    else
      j = j-1;
    end
    % transform to unit interval
    t = (x(k)-as(j))/(bs(j)-as(j));
    % evaluate the polynomial piece
    y(k) = interp_(xs,ws,ps(:,j),t);
  end
end

