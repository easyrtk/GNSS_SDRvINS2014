function dcmnb = eul2dcm(eulangle)
% conversion from Euler angles to direction cosine matrix
% input: eulrvec(1,2,3) = roll, pitch, yao angle /rad
% output  dcmnb = 3x3 direction cosine matrix from the navigation frame
%             to body frame
% reference:
% Titterton D H, Weston J L. Strapdown Inertial Navigation Technology[M], 2nd ed. 2004, p41.
% zsh 20140905

if nargin<1
    error('no input')
end
phi = eulangle(1);
theta = eulangle(2);
psi = eulangle(3);
cpsi = cos(psi);
spsi = sin(psi);
cthe = cos(theta);
sthe = sin(theta);
cphi = cos(phi);
sphi = sin(phi);

C1 = [cpsi  spsi 0;
      -spsi cpsi 0;
       0     0   1];
C2 = [cthe  0  -sthe;
      0   1     0 ;
      sthe  0   cthe];
C3 = [1   0    0;
      0  cphi sphi;
      0 -sphi cphi];  

dcmnb = C3 * C2 * C1;
