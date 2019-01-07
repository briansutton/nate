function varargout = plotpartition(varargin)
%PLOTPARTITION   Plot a partition of a real interval using vertical lines.
%
%   PLOTPARTITION([A B],N) plots the uniform partition of [A,B] consisting
%   of N subintervals.
%
%   PLOTPARTITION([ A1 A2 ... AN; B1 B2 ... BN ]) plots the partition
%   [AJ,BJ], J = 1,...,N.
%
%   Inputs (form 1):
%
%     AB = [A B]: endpoints of the entire interval
%
%     N: number of subintervals
%
%   Inputs (form 2):
%
%     ASBS = [ A1 A2 ... AN; B1 B2 ... BN ]: endpoints of the subintervals
%     [AJ,BJ], J = 1,...,N
%
%   Outputs:
%
%     H = PLOTPARTITION(...): graphics handle
%
%   Example:
%
%     f = @(x) cos(x); a = 0; b = 2*pi;
%     m = 2; n = 4;
%     xs = griduni([a b],m,n)
%     ps = sampleuni(f,[a b],m,n)
%     newfig;
%     plotfun(f,[a b]);
%     plotpartition([a b],n);
%     plotsample(xs,ps);
%
%   Copyright 2019 Brian Sutton

narginchk(1,2);
natecheck('plotpartition',varargin{:});
if nargin==1
  asbs = varargin{1};
  x = sort(asbs(:)');
else
  ab = varargin{1};
  n = varargin{2};
  a = min(ab);
  b = max(ab);
  x = linspace(a,b,n+1);
end
yl = ylim;
y = [ -logspace(20,-20,10)'; logspace(-20,20,10)' ];
h = plot(repmat(x,length(y),1),repmat(y,1,length(x))       ...
        ,'k','linewidth',0.5*get(0,'defaultlinelinewidth') ...
        ,'handlevisibility','off');
% Setting color with linespec and then changing color to desired RGB
% triplet prevents ColorOrderIndex from advancing in MATLAB R2017b.
set(h,'color',[0.5 0.5 0.5]);
xlim([x(1) x(end)]);
ylim(yl);
if nargout>0, varargout = { h }; end

