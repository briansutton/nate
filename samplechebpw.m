function ys = samplechebpw(f,asbs,m)
%SAMPLECHEBPW   Sample on a piecewise-Chebyshev grid.
%
%   YS = SAMPLECHEBPW(F,ASBS,M) evaluates a function at each node in a
%   piecewise-Chebyshev grid.
%
%   Inputs:
%
%     F: function
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     in a partition
%
%     M: degree
%
%   Outputs:
%
%     YS: matrix of function values
%
%   The piecewise-Chebyshev grid consisting of a Chebyshev grid of degree M
%   on [AJ,BJ], J = 1,...,N, is constructed. Then the function F is
%   evaluated at each node to produce YS.
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('samplechebpw',f,asbs,m);
xs = gridchebpw(asbs,m);
ys = arrayfun(f,xs);

