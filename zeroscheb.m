function z = zeroscheb(ps,ab,tol)
%ZEROSCHEB   Zeros of a polynomial interpolant on a Chebyshev grid.
%
%   Z = ZEROSCHEB(PS,[A B],TOL) computes the zeros of an interpolating
%   polynomial on a Chebyshev grid.
%
%   Inputs:
%
%     PS: Chebyshev sample of a polynomial on an interval [A,B]
%
%     AB = [A B]: endpoints of the interval
%
%     TOL: tolerance used in testing if a computed zero is close enough to
%     the given interval (default: TOL = 1e-6)
%
%   Outputs:
%
%     Z: zeros of the polynomial in the interval [A,B]
%
%   If the grid is of degree M and the polynomial is of degree at most M,
%   then the zeros are exact (in theory).
%
%   Example:
%
%     f = @(x) cos(x); a = 0; b = 10;
%     m = 30;
%     ps = samplecheb(f,[a b],m);
%     z = zeroscheb(ps,[a b])
%     newfig;
%     plotfun(f,[a b]);
%     plot(z,zeros(size(z)),'*');
%
%   Copyright 2019 Brian Sutton

narginchk(2,3);
if nargin<3, tol = 1e-6; end
natecheck('zeroscheb',ps,ab,tol);
a = ab(1);
b = ab(2);
m = size(ps,1)-1;
[xs,ws] = gridcheb([-1 1],m);
% compute generalized eigenvalues
e = zeros_(xs,ws,ps);
% select zeros in domain and discard other eigenvalues
e = e(abs(imag(e))<=tol&abs(real(e))<=1+tol);
z = sort((a+b)/2+(b-a)/2*e);

