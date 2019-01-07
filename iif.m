function z = iif(b,x,y)
%IIF   Inline if-then-else.
%
%   Z = IIF(B,X,Y) evaluates X or Y depending on whether B is true or
%   false.
%
%   Inputs:
%
%     B: a boolean condition
%
%     X: function to evaluate if B is true
%
%     Y: function to evaluate if B is false
%
%   Output:
%
%     Z: equals X() if B is true or Y() if B is false
%
%   Example:
%
%     x = -3;
%     iif(x>=0,@() x,@() -x)
%     y = 5;
%     iif(y>=0,@() y,@() -y)
%
%   Copyright 2019 Brian Sutton

if b
  z = x();
else
  z = y();
end

