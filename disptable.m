function disptable(varargin)
%DISPTABLE   Table display.
%
%   DISPTABLE(H1,V1,F1,H2,V2,F2,...) displays a table with column headings
%   and formatted entries.
%
%   Inputs:
%
%     H1: label for first column
%
%     V1: vector of entries for first column
%
%     F1: format string for first column
%
%     H2: label for second column
%
%     V2: vector of entries for second column
%
%     F2: format string for second column
%
%     ...
%
%   For more information on format strings, type HELP FPRINTF.
%
%   Example:
%
%     n = [0; 1; 2; 3; 4];
%     x = sqrt(n);
%     disptable('n',n,'%1d','sqrt(n)',x,'%7.4f');
%
%   Copyright 2019 Brian Sutton

narginchk(3,inf);
natecheck('disptable',varargin{:});

headings = varargin(1:3:end);
values = varargin(2:3:end);
formats = varargin(3:3:end);

m = max(cellfun(@numel,values));
n = nargin/3;

width = nan(1,n);
for j = 1:n
  width(j) = length(headings{j});
  for i = 1:length(values{j})
    s = sprintf([formats{j}],values{j}(i));
    if length(s)>width(j)
      width(j) = length(s);
    end
  end
end

for j = 1:n
  fprintf(['%' num2str(width(j)) 's'],headings{j});
  if j<n
    fprintf('   ');
  end
end
fprintf('\n');
for i = 1:m
  for j = 1:n
    if length(values{j})<i
      s = '';
    else
      s = sprintf([formats{j}],values{j}(i));
    end
    fprintf(['%-' num2str(width(j)) 's'],s);
    if j<n
      fprintf('   ');
    end
  end
  fprintf('\n');
end

