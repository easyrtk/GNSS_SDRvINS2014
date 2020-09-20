function plotDeep(begintime,endtime)
Rearth=6378388;
est_height1=[];
est_lon1=[];
est_lat1=[];
est_roll1=[];
est_pitch1=[];
est_yaw1=[];
load('C:\Sara\School\ASEN-6080StatOD2\Project1\inputs\SPANIMUreadings.mat')
if nargin==0
    endtime=1630;
    begintime=130;
end
timex=begintime+[0.1:0.1:endtime-130];
segno=(endtime-begintime)/100;

for i=1:segno
    fname=['deepresult',num2str(begintime+(i-1)*100),'~',num2str(begintime+i*100),'.mat'];
    load(fname);

    est_height1=[est_height1,est_height(2:end)];
    est_lon1=[est_lon1,est_lon(2:end)];
    est_lat1=[est_lat1,est_lat(2:end)];
    est_roll1=[est_roll1,est_roll_KF(2:end)];
    est_pitch1=[est_pitch1,est_pitch_KF(2:end)];
    est_yaw1=[est_yaw1,est_yaw_KF(2:end)];
end
npts=(endtime-begintime)*100+1;

figure
plot((INSlon(ind1:10:ind1+npts-1)/180*pi+1.8375)*Rearth*cos(0.6984),(INSlat(ind1:10:ind1+npts-1)/180*pi-0.6984)*Rearth)
hold on
plot((est_lon1(1:100:end)+1.8375)*Rearth*cos(0.6984),(est_lat1(1:100:end)-0.6984)*Rearth,'r')
title('trajectory')
hold off
errlat=(est_lat1(1:10:end)'-INSlat(ind1:ind1+npts-2)/180*pi)*Rearth;
errlon=(est_lon1(1:10:end)'-INSlon(ind1:ind1+npts-2)/180*pi)*Rearth*cos(0.6984);
errhei=est_height1(1:10:end)'-INShei(ind1:ind1+npts-2);

figure
plot(timex,sqrt(errlat(1:10:end).^2+errlon(1:10:end).^2));
title('horizontal error')

% figure
% plot(timex,INShei(ind1:10:ind1+npts-2))
% hold on
% plot(timex,est_height1(1:100:end),'r')
% title('height')
% hold off

figure
plot(timex,errhei(1:10:end));
title('vertical error');

figure
plot(timex,sqrt(errlat(1:10:end).^2+errlon(1:10:end).^2+errhei(1:10:end).^2))
title('3D error')
errheiVDI=errhei;
errlatVDI=errlat;
errlonVDI=errlon;
est_heightVDI=est_height1;
est_latVDI=est_lat1;
est_lonVDI=est_lon1;
timexVDI=timex;
save ('deepsolutions.mat','errheiVDI','errlatVDI','errlonVDI','est_heightVDI','est_latVDI','est_lonVDI','timexVDI')
fprintf('rms error: lat:%.2f, lon:%.2f, h:%.2f\n',rms(errlatVDI),rms(errlonVDI),rms(errheiVDI))
errroll=(est_roll1(1:10:end)'/pi*180-INSroll(ind1-111:ind1-111+npts-2));
errpitch=(est_pitch1(1:10:end)'/pi*180-INSpitch(ind1-111:ind1-111+npts-2));
erryaw0=est_yaw1(1:10:end)'/pi*180-INShead(ind1-111:ind1-111+npts-2);
modyaw=zeros(length(errroll),1);
modyaw(find(erryaw0<-360))=-360;
modyaw(((erryaw0>=-360)&(erryaw0<=-180))|(erryaw0>=0))=360;
modyaw((erryaw0>=-180)&(erryaw0<0))=-360;
erryaw=mod(erryaw0,modyaw);

figure
subplot(3,1,1)
plot(timex,errroll(1:10:end));
title('roll error');
subplot(3,1,2)
plot(timex,errpitch(1:10:end));
title('pitch error');
subplot(3,1,3)
plot(timex,erryaw(1:10:end));
title('yaw error');