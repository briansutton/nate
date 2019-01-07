function [M,x] = infnorm(f,ab,n)
%INFNORM   Infinity norm of a function. (This is also known as the supremum
%norm.)
%
%   [M,X] = INFNORM(F,[A B],N) computes the infinity norm of the function F
%   over the interval [A,B].
%
%   Inputs:
%
%     F: a real-valued function of a real variable
%
%     AB = [A B]: interval endpoints
%
%     N: number of sample points to start
%
%   Outputs:
%
%     M: infinity norm of F on the interval [A,B]
%
%     X: location of a maximum of ABS(F(X))
%
%   The search for the maximum absolute value begins with N points. Then
%   the approximate maximum is refined using a local search. Larger values
%   of N generally lower the chance of missing a spike in the function.
%
%   The search spends a considerable portion of its resources near X = A
%   and X = B because singularities lie at the endpoints in some common
%   scenarios. A spike in the middle of the interval may easily be missed
%   when N is small.
%
%   Example:
%
%     f = @(x) sin(7*x)+cos(11*x); a = -1; b = 1;
%     M = infnorm(f,[a b],100)
%     newfig;
%     plotfun(@(x) f(x),[a b]);
%     plotfun(@(x) M,[a b]);
%     plotfun(@(x) -M,[a b]);
%
%   Copyright 2019 Brian Sutton

narginchk(3,3);
natecheck('infnorm',f,ab,n);
c = ab(1);
d = ab(2);
x = mean([c d])-(d-c)/2*sin(pi/2*cos(pi*randpoints_(n)));
y = arrayfun(f,x);
[~,i] = max(abs(y));
if i==1, c = c; else c = x(i-1); end
if i==length(x), d = d; else d = x(i+1); end
if c<x(i)&&~isfinite(f(c))
  cc = c;
  dd = x(i);
  for j = 1:10
    if isfinite(f((cc+dd)/2))
      dd = (cc+dd)/2;
    else
      cc = (cc+dd)/2;
    end
  end
  c = dd;
end
if d>x(i)&&~isfinite(f(d))
  cc = x(i);
  dd = d;
  for j = 1:10
    if isfinite(f((cc+dd)/2))
      cc = (cc+dd)/2;
    else
      dd = (cc+dd)/2;
    end
  end
  d = cc;
end
[M,x] = golden_(@(x) abs(f(x)),c,d);

