function ylog
%YLOG   Set the vertical axis to logarithmic.
%
%   YLOG sets the scale on the vertical axis of the current plot to be
%   logarithmic.
%
%   Copyright 2019 Brian Sutton

ax = gca;
ax.YScale = 'log';
hold on;

