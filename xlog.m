function xlog
%XLOG   Set the horizontal axis to logarithmic.
%
%   XLOG sets the scale on the horizontal axis of the current plot to be
%   logarithmic.
%
%   Copyright 2019 Brian Sutton

ax = gca;
ax.XScale = 'log';
hold on;

