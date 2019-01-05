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
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('resamplematrixcheb',l,m);
[xs,ws] = gridcheb([-1 1],m);
xs_ = gridcheb([-1 1],l);
E = evalmatrix_(xs_,xs,ws);

