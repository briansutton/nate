function p = interpchebpw(ps,asbs)
%INTERPCHEBPW   Interpolation on a piecewise-Chebyshev grid.
%
%   P = INTERPCHEBPW(PS,ASBS) constructs a piecewise-polynomial interpolant
%   P with given values on a piecewise-Chebyshev grid. The subintervals of
%   the piecewise grid do not have to be of equal width.
%
%   Inputs:
%
%     PS: an (M+1)-by-N matrix containing values of the interpolant on the
%     piecewise-Chebyshev grid
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     [AJ,BJ] in the piecewise-defined grid
%
%   Outputs:
%
%     P: the interpolating polynomial
%
%   Let XS be the piecewise-Chebyshev grid consisting of a degree-M
%   Chebyshev grid on each subinterval. Then P is constructed so that
%   P(XS(I,J)) = PS(I,J). In addition, P is a polynomial of degree at most
%   M on each subinterval [AJ,BJ].
%
%   Example:
%
%     f = @(x) exp(x)*sin(3*x); a = -pi; b = pi;
%     asbs = [ -pi 0 pi/2; 0 pi/2 pi ];
%     m = 4;
%     ps = samplechebpw(f,asbs,m);
%     p = interpchebpw(ps,asbs);
%     newfig;
%     plotfun(f,[a b]);
%     plotpartition(asbs);
%     plotsample(gridchebpw(asbs,m),ps);
%     plotfun(p,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('interpchebpw',ps,asbs);
m = size(ps,1)-1;
n = size(ps,2);
% arrange subintervals into ascending order
for j = 1:n
  if asbs(1,j)>asbs(2,j)
    asbs(:,j) = flipud(asbs(:,j));
    ps(:,j) = flipud(ps(:,j));
  end
end
[~,sigma] = sort(asbs(1,:));
as = asbs(1,sigma);
bs = asbs(2,sigma);
ps = ps(:,sigma);
% construct the interpolant
[xs,ws] = gridcheb([0 1],m);
p = @(x) interpchebpw_(xs,ws,ps,as,bs,x);

