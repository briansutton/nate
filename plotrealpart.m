function varargout = plotrealpart(f,lims)
%PLOTREALPART   Surface plot of the real part of a complex function.
%
%   PLOTREALPART(F,LIMS) plots REAL(F(X+I*Y)) as a surface in 3-D.
%
%   Inputs:
%
%     F: a function from complex numbers to complex numbers
%
%     LIMS = [XMIN XMAX YMIN YMAX ZMIN ZMAX]: axes limits
%
%   Outputs:
%
%     H = PLOTIMAGPART(...): graphics handle
%
%   Example:
%
%     f = @(z) sqrt(z);
%     xmin = -1; xmax = 1;
%     ymin = -1; ymax = 1;
%     vmin = -1.5; vmax = 1.5;
%     newfig;
%     plotrealpart(f,[xmin xmax ymin ymax vmin vmax]);
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('plotrealpart',f,lims);
s = warning('off','MATLAB:fplot:NotVectorized');
legend hide;
h = fsurf(@(x,y) real(f(x+1i*y)),lims(1:4));
zlim(lims(5:6));
caxis(lims(5:6));
view(3);
warning(s);
if nargout>0, varargout = { h }; end

