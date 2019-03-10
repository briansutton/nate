function varargout = plotfun(f,ab,varargin)
%PLOTFUN   Plot a function of a continuous variable.
%
%   PLOTFUN(F,[A B],...) plots a function F(X) on an interval [A,B].
%
%   Inputs:
%
%     F: function to plot
%
%     AB = [A B]: endpoints of interval
%
%     ...: additional arguments for PLOT
%
%   Outputs:
%
%     H = PLOTFUN(...): graphics handle to the line object
%
%   Additional arguments beyond AB control formatting, the legend entry,
%   etc. See HELP PLOT.
%
%   The resolution of the plot is reasonably high relative to the width of
%   the interval [A,B]. If a zoomed-in graph on some smaller interval [C,D]
%   is desired, then it may be necessary to call PLOTFUN again rather than
%   simply resetting the zoom window using XLIM or AXIS.
%
%   Example:
%
%     newfig;
%     plotfun(@(x) exp(x),[-2 2],'displayname','exp');
%     xlabel('x');
%
%   Copyright 2019 Brian Sutton

narginchk(2,inf);
natecheck('plotfun',f,ab);
a = ab(1);
b = ab(2);
n = 250;
x = a+(b-a)*randpoints_(n);
y = arrayfun(f,x);
xx = nan(n-1,1);
yy = nan(n-1,1);
m = 0;
for k = 1:n-2
  if y(k+1)>y(k)&&y(k+1)>y(k+2)
    m = m+1;
    [yy(m),xx(m)] = golden_(f,x(k),x(k+2));
  elseif y(k+1)<y(k)&&y(k+1)<y(k+2)
    m = m+1;
    [yy(m),xx(m)] = golden_(@(x) -f(x),x(k),x(k+2));
    yy(m) = -yy(m);
  elseif isfinite(y(k))&&~isfinite(y(k+1))
    m = m+1;
    c = x(k);
    d = x(k+1);
    for i = 1:10
      if isfinite(f((c+d)/2))
        c = (c+d)/2;
      else
        d = (c+d)/2;
      end
    end
    xx(m) = c;
    yy(m) = f(c);
  elseif ~isfinite(y(k+1))&&isfinite(y(k+2))
    m = m+1;
    c = x(k+1);
    d = x(k+2);
    for i = 1:10
      if isfinite(f((c+d)/2))
        d = (c+d)/2;
      else
        c = (c+d)/2;
      end
    end
    xx(m) = d;
    yy(m) = f(d);
  end
end
xxx = linspace(a,b,1000)';
x = [x; xx(1:m); xxx];
[x,I] = sort(x);
y = [y; yy(1:m); arrayfun(f,xxx)];
y = y(I);
y(y==0) = realmin;
h = plot(x,y,varargin{:});
xlim(sort(ab));
if nargout>0, varargout = { h }; end

