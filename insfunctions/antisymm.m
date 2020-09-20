function mat = antisymm(vect)
% convert a 3x1 vector into anti-symmetric matrix. 
% input:    vect(1,2,3)=x,y,z
% output:
%       mat = [ 0 -z  y;
%             z 0  -x;
%            -y x   0]
% zsh 20140901

x=vect(1);
y=vect(2);
z=vect(3);
mat=[0 -z  y
     z  0 -x
     -y x  0];
  
  
