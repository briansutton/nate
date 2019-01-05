function [as,bs] = partition_(ab,n)
%PARTITION_   Divide an interval into subintervals of equal width.
%
%   [AS,BS] = PARTITION_([A B],N) partitions the interval [A,B] into N
%   contiguous subintervals of equal width.
%
%   Inputs:
%
%     AB = [A B]: endpoints for the original interval
%
%     N: number of subintervals
%
%   Outputs:
%
%     AS: vector of left endpoints
%
%     BS: vector of right endpoints
%
%   The subintervals are [AS(J),BS(J)], J = 1,...,N.
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('partition_',ab,n);
a = ab(1);
b = ab(2);
endpoints = linspace(a,b,n+1);
as = endpoints(1:end-1);
bs = endpoints(2:end);

