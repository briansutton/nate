function y = polynomialeval(c,x)
%POLYNOMIALEVAL   Polynomial evaluation.
%
%   POLYNOMIALEVAL evaluates a polynomial of the form
%
%     P(X) = C0 + C1*X + C2*X^2 + ... + CN*X^N
%
%   at a given value of X.
%
%   Inputs:
%
%     C: vector of coefficients [C0 ... CN]
%
%     X: location where the polynomial is evaluated
%
%   Outputs:
%
%     Y: value of the polynomial at the given X
%
%   Example:
%
%     polynomialeval([2 -1 4],0.5)
%     2-0.5+4*0.5^2
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('polynomialeval',c,x);
if isempty(c), y = 0; return; end
n = length(c)-1;
y = c(n+1);
for j = n:-1:1
  y = c(j)+x*y;
end

