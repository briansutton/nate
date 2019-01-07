function [z,x] = newton(f,fprime,x0,kmax,tol)
%NEWTON   Newton iteration for locating a zero of a function.
%
%   [Z,X] = NEWTON(F,FPRIME,X0,KMAX,TOL) solves F(X) = 0 using Newton's
%   method.
%
%   Inputs:
%
%     F: a function
%
%     FPRIME: the derivative of F
%
%     X0: a starting guess for the zero
%
%     KMAX: maximum number of Newton steps (default: KMAX = 100)
%
%     TOL: termination tolerance (default: TOL = 1e-13)
%
%   Outputs:
%
%     Z: computed solution
%
%     X: sequence of Newton iterates culminating in Z
%
%   If Newton's method fails to converge, then Z is assigned NaN.
%
%   Example:
%
%     f = @(x) log(x)-1;
%     fprime = @(x) 1/x;
%     x0 = 2.8;
%     [z,x] = newton(f,fprime,x0)
%     newfig;
%     plotfun(f,[0.1 4]);
%     axis([0 4 -3 3]);
%     grid on;
%     plot(z,0,'*');
%
%   Copyright 2019 Brian Sutton

narginchk(3,5);
if nargin<4, kmax = 100; end
if nargin<5, tol = 1e-13; end
natecheck('newton',f,fprime,x0,kmax,tol);
x = nan(kmax+1,1);
% initialize
x(1) = x0;
% iterate
for k = 1:kmax
  fx = f(x(k));
  % terminate if solution found
  if fx==0
    x = x(1:k); z = x(end); return;
  end
  % take a Newton step
  fpx = fprime(x(k));
  dx = -fx/fpx;
  x(k+1) = x(k)+dx;
  % signal failure on NaN or infinity
  if ~isfinite(x(k+1))
    x = x(1:k+1); z = NaN; return;
  end
  % terminate on numerical convergence
  if abs(dx)<=tol*max(abs(x(k+1)),1)
    x = x(1:k+1); z = x(end); return;
  end
end
% signal failure to converge
z = NaN;

