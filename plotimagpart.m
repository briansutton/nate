function varargout = plotimagpart(f,lims)
%PLOTIMAGPART   Surface plot of the imaginary part of a complex function.
%
%   PLOTIMAGPART(F,LIMS) plots IMAG(F(X+I*Y)) as a surface in 3-D.
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
%     plotimagpart(f,[xmin xmax ymin ymax vmin vmax]);
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('plotimagpart',f,lims);
s = warning('off','MATLAB:fplot:NotVectorized');
legend hide;
h = fsurf(@(x,y) imag(f(x+1i*y)),lims(1:4));
zlim(lims(5:6));
caxis(lims(5:6));
view(3);
warning(s);
if nargout>0, varargout = { h }; end

