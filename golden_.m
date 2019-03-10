function [y,x] = golden_(f,a,b)
%GOLDEN_   Locate a local maximum using golden-section search.
%
%   [Y,X] = GOLDEN_(F,A,B) finds a local maximum of F in [A,B].
%
%   Inputs:
%
%     F: a continuous real-valued function of a real variable
%
%     A: left endpoint of search interval
%
%     B: right endpoint of search interval
%
%   Outputs:
%
%     Y: local maximum value Y = F(X)
%
%     X: location of a local maximum
%
%   This is an internal routine called by INFNORM and PLOTFUN.
%
%   Copyright 2019 Brian Sutton

c = min(a,b);
d = max(a,b);
phi = 1.618033988749895;
n = ceil((log(d-c)-log(max(eps(a),eps(b))))/log(phi));
delta = (d-c)/phi;
x1 = d-delta;
x2 = c+delta;
f1 = f(x1);
f2 = f(x2);
for k = 1:n
  if f1>f2
    d = x2;
    delta = (d-c)/phi;
    x2 = x1;
    x1 = d-delta;
    f2 = f1;
    f1 = f(x1);
  else
    c = x1;
    delta = (d-c)/phi;
    x1 = x2;
    x2 = c+delta;
    f1 = f2;
    f2 = f(x2);
  end
end
x = [ a c x1 x2 d b ];
y = [ f(a) f(c) f1 f2 f(d) f(b) ];
[y,i] = max(y);
x = x(i);
