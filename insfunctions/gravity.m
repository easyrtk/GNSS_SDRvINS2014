function g = gravity(lat,height)
%calculate gravity 
%       
%   input
%       lat = latitude /rad
%       height = height above the ellipsoid /m
%
%  output
%       g: gravity / m/(s*s)
%
% reference
% Titterton D H, Weston J L. Strapdown Inertial Navigation Technology[M], 2nd ed. 2004 p57.
% zsh 20140901

radiusa = localradius(lat);
g0 = 9.780318*( 1 + 5.3024e-3*(sin(lat))^2 - 5.9e-6*(sin(2*lat))^2 );
g = g0/(1 + height/radiusa)^2;
if height < -radiusa
  error('Input height out of boundary')
end
   