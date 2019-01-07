function ylog
%YLOG   Set the vertical axis to logarithmic.
%
%   YLOG sets the scale on the vertical axis of the current plot to be
%   logarithmic.
%
%   Example:
%
%     newfig;
%     plotfun(@(x) exp(-x),[0 10]);
%     ylog;
%
%   Copyright 2019 Brian Sutton

ax = gca;
ax.YScale = 'log';
hold on;

