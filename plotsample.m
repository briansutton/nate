function varargout = plotsample(xs,ys,varargin)
%PLOTSAMPLE   Plot a sequence of points.
%
%   PLOTSAMPLE(XS,YS,...) plots points (XS(J),YS(J)), J = 1,...,N.
%
%   Inputs:
%
%     XS: horizontal coordinates
%
%     YS: vertical coordinates
%
%     ...: additional arguments for PLOT
%
%   Outputs:
%
%     H = PLOTSAMPLE(...): graphics handle
%
%   Additional arguments beyond YS control formatting, the legend entry,
%   etc. See HELP PLOT.
%
%   Points are plotted in a single series, regardless of whether they are
%   provided as row or column vectors or matrices.
%
%   Example:
%
%     f = @(x) exp(x);
%     xs = [-2; -1; 0; 1; 2];
%     ys = arrayfun(f,xs);
%     newfig;
%     plotfun(f,[-2 2]);
%     plotsample(xs,ys);
%
%   Copyright 2019 Brian Sutton

narginchk(2,inf);
natecheck('plotsample',xs,ys);
if nargin<3
  h = plot(xs(:),ys(:),'k.' ...
          ,'markersize',3*get(0,'defaultlinemarkersize'));
  % Setting color with linespec and then changing color to desired RGB
  % triplet prevents ColorOrderIndex from advancing in MATLAB R2017b.
  h.MarkerEdgeColor = [0.3 0.3 0.3];
else
  h = plot(xs(:),ys(:),varargin{:});
end
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
if nargout>0, varargout = { h }; end

