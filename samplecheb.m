function ys = samplecheb(f,ab,m)
%SAMPLECHEB   Sample on a Chebyshev grid.
%
%   YS = SAMPLECHEB(F,[A B],M) evaluates a function at each node in a
%   Chebyshev grid.
%
%   Inputs:
%
%     F: function
%
%     AB = [A B]: endpoints of an interval
%
%     M: degree
%
%   Outputs:
%
%     YS: vector of function values
%
%   YS(I) = F(XS(I)), in which XS is the Chebyshev grid of degree M on
%   [A,B].
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('samplecheb',f,ab,m);
xs = gridcheb(ab,m);
ys = arrayfun(f,xs);

