function p = interpgen(xs,ps)
%INTERPGEN   Interpolation on a given grid.
%
%   P = INTERPGEN(XS,PS) constructs an interpolating polynomial P whose
%   graph passes through given points.
%
%   Inputs:
%
%     XS: horizontal coordinates of interpolation points
%
%     PS: vertical coordinates of interpolation points
%
%   Outputs:
%
%     P: the interpolating polynomial
%
%   No two entries of XS may be equal. If XS is (M+1)-by-1, then the
%   entries form a grid of degree M, and there is a unique interpolating
%   polynomial P of degree at most M.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('interpgen',xs,ps);
ws = baryweights(xs);
p = @(x) interp_(xs,ws,ps,x);

