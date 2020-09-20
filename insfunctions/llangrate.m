function omega_el_L = llangrate(lat,ve,vn,height)
%  calculate angular rate of the local-level frame, e.g. ENU 
%
%  input
%       lat = latitude/rad
%       ve = east velocity / m/s
%       vn = north velocity / m/s
%       height = height above the reference ellipsoid /m
%
%  output
%       omega_el_l: local level frame angular rate vector relative to ecef, expressed in local-level frame
%           omega_el_l(1) = local-level frame x axis angular rate /(rad/s)
%           omega_el_l(2) = local-level frame y axis angular rate /(rad/s)
%           omega_el_l(3) = local-level frame z axis angular rate /(rad/s)
%
%  reference: 董绪荣, 张守信, 华仲春, GPS/INS组合导航定位及其应用 p64, ISBN: 7-81024-472-8
%  Titterton D H, Weston J L. Strapdown Inertial Navigation Technology[M], 2nd ed. 2004, p49.
%  zsh 20140902

[~,rm,rp]=localradius(lat);
omega_el_L(1,1)=-vn/(rm+height);
omega_el_L(2,1)=ve/(rp+height);
omega_el_L(3,1)=ve*tan(lat)/(rp+height);% could be set to zero
