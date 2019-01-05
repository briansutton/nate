function z = zerosgen(xs,ps)
%ZEROSGEN   Zeros of a polynomial interpolant.
%
%   Z = ZEROSGEN(XS,PS) computes the zeros of an interpolating polynomial.
%
%   Inputs:
%
%     XS: horizontal coordinates of interpolation points
%
%     PS: vertical coordinates of interpolation points
%
%   Outputs:
%
%     Z: zeros of the polynomial interpolant
%
%   If XS is a grid of degree M, then the polynomial interpolant is the
%   unique polynomial of degree at most M whose graph intersects the given
%   points, and the computed zeros are the exact zeros of the polynomial
%   (in theory).
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('zerosgen',xs,ps);
ws = baryweights(xs);
% compute generalized eigenvalues
z = zeros_(xs,ws,ps);
% discard infinite eigenvalues, leaving polynomial zeros
z = z(isfinite(z));

