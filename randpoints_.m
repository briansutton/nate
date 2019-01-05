function x = randpoints_(n)
%RANDPOINTS_   Quasi-random numbers in [0,1].
%
%   RANDPOINTS_ generates random numbers that are almost uniformly spaced
%   in [0,1].
%
%   Inputs:
%
%     N: number of points
%
%   Outputs:
%
%     X: vector of N random points in [0,1]
%
%   This is an internal routine called by INFNORM and PLOTFUN.
%
%   Copyright 2019 Brian Sutton

persistent points
natecheck('randpoints_',n);
if n<=0, x = []; return; end
if n>10000
  x = unique([0; linspace(0,1-1/n,n-2)'+1/n*rand(n-2,1); 1]);
  return;
end
if isempty(points)
  points = rand(RandStream('mt19937ar','Seed',0),10000,1);
end
x = unique([0; linspace(0,1-1/n,n-2)'+1/n*points(1:n-2,1); 1]);

