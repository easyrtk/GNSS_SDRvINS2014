function [radiusa,rm,rp] = localradius(lat)
% calculate meridian and prime radius 
%   
%  input
%       lat: geodetic latitude /rad
%
%  output
%       radiusa: mean curvature radius at input latitude /m
%       rm: meridian radius /m
%       rp: prime radius /m
% reference:
% 董绪荣, 张守信, 华仲春, GPS/INS组合导航定位及其应用 pp63-64 ISBN: 7-81024-472-8
% Titterton D H, Weston J L. Strapdown Inertial Navigation Technology[M], 2nd ed. 2004, p49.
% zsh 20140902


a = 6378137;    % major semi axis/m
e2 = 0.0066943799901378; % square of eccentricity ecc = 0.0818191908426
tmp = 1-e2*(sin(lat))^2;
rm = (a*(1-e2))/( (tmp)^(3/2) );
rp = a/( sqrt(tmp) );
radiusa=sqrt(rm*rp);