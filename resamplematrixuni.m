function E = resamplematrixuni(l,m)
%RESAMPLEMATRIXUNI   Matrix for resampling from one uniform grid to
%another.
%
%   E = RESAMPLEMATRIXUNI(L,M) constructs a resampling matrix for uniform
%   grids.
%
%   Inputs:
%
%     L: degree of output grid
%
%     M: degree of input grid
%
%   Outputs:
%
%     E: resampling matrix
%
%   Suppose P is a polynomial of degree at most M and PS is a sample of P
%   on a uniform grid of degree M. Then E*PS is a sample of P on the
%   uniform grid of degree L over the same interval.
%
%   Example:
%
%     l = 4;
%     m = 3;
%     E = resamplematrixuni(l,m);
%     p = @(x) 2*x^3-x^2+4*x-3; a = -1; b = 1;
%     ps = sampleuni(p,[a b],m);
%     newfig;
%     plotfun(p,[a b]);
%     plotsample(griduni([a b],m),ps);
%     plotsample(griduni([a b],l),E*ps,'o');
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('resamplematrixuni',l,m);
[xs,ws] = griduni([0 1],m,1);
xs_ = griduni([0 1],l,1);
E = evalmatrix_(xs_,xs,ws);

