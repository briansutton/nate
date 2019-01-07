function varargout = plotfun3(fx,fy,fz,ab,varargin)
%PLOTFUN3   Plot a parametric function as a curve in 3-D.
%
%   PLOTFUN3(FX,FY,FZ,[A B],...) plots a parameterized path
%   (FX(T),FY(T),FZ(T)) from T = A to T = B.
%
%   Inputs:
%
%     FX: X-coordinate function
%
%     FY: Y-coordinate function
%
%     FZ: Z-coordinate function
%
%     AB = [A B]: limits for the parameter T
%
%     ...: additional arguments for FPLOT3
%
%   Outputs:
%
%     H = PLOTFUN3(...): graphics handle
%
%   Additional arguments beyond AB control formatting, the legend entry,
%   etc. See HELP FPLOT3.
%
%   Example:
%
%     newfig;
%     plotfun3(@(t) cos(t),@(t) sin(t),@(t) t,[0 8*pi]);
%
%   Copyright 2019 Brian Sutton

narginchk(4,inf);
natecheck('plotfun3',fx,fy,fz,ab);
s = warning('off','MATLAB:fplot:NotVectorized');
h = fplot3(fx,fy,fz,sort(ab),varargin{:});
[az,el] = view;
if el==90, view(3); end
warning(s);
if nargout>0, varargout = { h }; end

