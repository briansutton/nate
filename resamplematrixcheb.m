function E = resamplematrixcheb(l,m)
%RESAMPLEMATRIXCHEB   Matrix for resampling from one Chebyshev grid to
%another.
%
%   E = RESAMPLEMATRIXCHEB(L,M) constructs a Chebyshev resampling matrix.
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
%   on a Chebyshev grid of degree M. Then E*PS is a sample of P on the
%   Chebyshev grid of degree L over the same interval.
%
%   Example:
%
%     l = 4;
%     m = 3;
%     E = resamplematrixcheb(l,m);
%     p = @(x) 2*x^3-x^2+4*x-3; a = -1; b = 1;
%     ps = samplecheb(p,[a b],m);
%     newfig;
%     plotfun(p,[a b]);
%     plotsample(gridcheb([a b],m),ps);
%     plotsample(gridcheb([a b],l),E*ps,'o');
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('resamplematrixcheb',l,m);
[xs,ws] = gridcheb([-1 1],m);
xs_ = gridcheb([-1 1],l);
E = evalmatrix_(xs_,xs,ws);

