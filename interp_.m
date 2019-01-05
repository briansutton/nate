function px = interp_(xs,ws,ys,x)
%INTERP_   Interpolation using the barycentric form.
%
%   PX = INTERP_(XS,WS,YS,X) interpolates the points (XS,YS) with a
%   polynomial P and then evaluates P(X).
%
%   Inputs:
%
%     XS: horizontal coordinates of interpolation points
%
%     WS: barycentric weights associated with the grid XS
%
%     YS: vertical coordinates of interpolation points
%
%     X: location where the interpolating polynomial should be evaluated
%
%   Outputs:
%
%     PX: the value P(X) of the interpolating polynomial at the desired
%     location
%
%   This is an internal routine called by INTERPGEN, INTERPUNI, and
%   INTERPCHEB.
%
%   Copyright 2019 Brian Sutton

px = nan(size(x));
for k = 1:numel(x)
  % evaluate barycentric formula
  zs = ws./(x(k)-xs);
  px(k) = sum(ys.*zs)/sum(zs);
  % fix if at node
  if isnan(px(k))
    i = find(~isfinite(zs),1);
    px(k) = ys(i);
  end
end

