function eulangle = dcm2eul(DCMbn)
%  calculate Euler angles from direction cosine matrix
%       
%  input
%     DCMbn = 3x3 direction cosine matrix transforming from body to
%     navigation frame
%
%  output
%     eulangle(1) = roll /rad
%     eulangle(2) = pitch /rad
%     eulangle(3) = yaw /rad
%
%  note: the case that pitch (theta) is close to +/-90 deg is not
%  considered here
%
%  reference:
%  Titterton D H, Weston J L. Strapdown Inertial Navigation Technology[M], 2nd ed. 2004, p46.
% zsh 20140905
phi=atan(DCMbn(3,2)/DCMbn(3,3));
theta=asin(-DCMbn(3,1));
psi=atan2(DCMbn(2,1),DCMbn(1,1));
eulangle=[phi,theta,psi]';
