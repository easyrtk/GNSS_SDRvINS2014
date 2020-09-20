function DCM=calDCM(ang_vect)
% input
%       ang_vect: angles/rad
% output
%       DCM: direction cosine matrix
% reference: Titterton D H, Weston J L. Strapdown Inertial Navigation Technology[M], 2nd ed. 2004, p312.
% zsh 20140903
as = antisymm(ang_vect);
mag = norm(ang_vect);
if mag == 0
  DCM = eye(3);
else
  DCM = eye(3)+(sin(mag)/mag)*as+((1-cos(mag))/mag^2)*as*as;
end