function z = zeros_(xs,ws,ps)
%ZEROS_   Computational kernel for computing zeros of a polynomial.
%
%   Z = ZEROS_(XS,WS,PS) computes all zeros Z of a polynomial.
%
%   Inputs:
%
%     XS: interpolation grid
%
%     WS: barycentric weights for XS
%
%     PS: sample of a function on XS
%
%   Outputs:
%
%     Z: zeros of the polynomial
%
%   If XS is of degree M and the sampled function is a polynomial of degree
%   at most M, then the zeros are computed exactly (in theory).
%
%   This is an internal routine called by ZEROSGEN, ZEROSUNI, and
%   ZEROSCHEB.
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('zeros_',xs,ws,ps);
if max(abs(ps))==0, z = mean([ min(xs) max(xs) ]); return; end
m = length(ps)-1;
A = [ 0 -ps'; ws diag(xs) ];
B = blkdiag(0,eye(m+1));
z = eig(A,B);

