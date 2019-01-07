function varargout = newfig(m,n)
%NEWFIG   New figure.
%
%   NEWFIG clears the current figure and optionally partitions the figure
%   window into an M-by-N grid of axes.
%
%   Inputs:
%
%     M: number of rows of plots (default: M = 1)
%
%     N: number of columns of plots (default: N = 1)
%
%   Outputs:
%
%     AX: an M-by-N grid of blank axes
%
%   The current figure is cleared to create the new figure. Execute
%   SUBPLOT(AX(I,J)) to activate the axes in the I,J position.
%
%   Example:
%
%     ax = newfig(1,3);
%     plotfun(@(x) 1,[-2 2]);
%     ylim([-4 4]);
%     subplot(ax(1,2));
%     plotfun(@(x) x,[-2 2]);
%     ylim([-4 4]);
%     subplot(ax(1,3));
%     plotfun(@(x) x^2,[-2 2]);
%     ylim([-4 4]);
%
%   Copyright 2019 Brian Sutton

narginchk(0,2);
if nargin<1
  m = 1;
  n = 1;
elseif nargin<2
  n = 1;
end
natecheck('newfig',m,n);
clf;
ax = nan(m,n);
for i = 1:m
  for j = 1:n
    h = subplot(m,n,(i-1)*n+j);
    ax(i,j) = h;
    h.ColorOrder = [ 0     0.447 0.698 ;
                     0.835 0.369 0     ;
                     0     0.620 0.451 ;
                     0.8   0.475 0.655 ];
    set(h,'defaultLineLineWidth',1.6);
    set(h,'defaultParameterizedFunctionLineLineWidth',1.6);
    h.FontSize = 14;
    box on;
    hold on;
    legend show;
  end
end
subplot(ax(1,1));
if nargout>0, varargout = { ax }; end

