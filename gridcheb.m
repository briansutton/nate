function [xs,ws] = gridcheb(ab,m)
%GRIDCHEB   Chebyshev grid.
%
%   [XS,WS] = GRIDCHEB([A B],M) computes the degree-M Chebyshev grid XS on
%   [A,B] and the corresponding barycentric weights WS.
%
%   Inputs:
%
%     AB = [A B]: endpoints of an interval
%
%     M: degree of the grid
%
%   Outputs:
%
%     XS: the Chebyshev grid of degree M over [A,B]
%
%     WS: barycentric weights for XS
%
%   Example:
%
%     a = 0; b = 2;
%     m = 6;
%     [xs,ws] = gridcheb([a b],m)
%
%   Copyright 2019 Brian Sutton

narginchk(2,2);
natecheck('gridcheb',ab,m);
a = ab(1);
b = ab(2);
xs = [ a; (a+b)/2-(b-a)/2*cos((1:m-1)'*pi/m); b ];
ws = [ 1/2; (-1).^(1:m-1)'; (1/2)*(-1)^m ];

