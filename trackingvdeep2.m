function [trackResults, channel]= trackingvdeep2(fid, channel,trackRes,navSolutions,eph,activeChnList,svTimeTable, settings)
% Performs code and carrier tracking  and navigation based on INS/GPS deep integration
%
%[trackResults, channel] = tracking(fid, channel,trackRes,navSolutions,eph,activeChnList,svTimeTable, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record for I
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       trackRes        -tracking results using scalar loop to initialize
%       navSolusions    -navigation solutions from scalar loop to
%                       initialize
%       eph             - ephemerides
%       activeChnList   - a list of active satellites in the dataset
%       svTimeTable     - satellite time to find transmit time of a sample
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
% modified to integrate ins and gps deeply -zsh 2014.02

%% Initialize result structure ============================================
tic
fprintf('start time:%s\n',datestr(now));
kpt=1e-3;%kalman filter process update time
kmt=1e-3;%kalman filter measurement update time/s
pdi=1e-3;%integration time /s
r2d=180/pi;
d2r=pi/180;
StartTime=130000;%start time in dataset /ms
begintime=StartTime/1000;% start time for plots in the end /s
tracklengthall=1500000/(pdi*1000); %TOTAL TRACKING LENGTH in ms
tracklength=100000/(pdi*1000);% vector tracking length FOR ONE SEGMENT ms
endtime=begintime+tracklengthall/1000;% end time for plots in the end /s
sectno=fix(tracklengthall/tracklength);

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, tracklength);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, tracklength);
trackResults.remCodePhase = zeros(1, tracklength);
% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, tracklength);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, tracklength);
trackResults.I_E            = zeros(1, tracklength);
trackResults.I_L            = zeros(1, tracklength);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, tracklength);
trackResults.Q_P            = zeros(1, tracklength);
trackResults.Q_L            = zeros(1, tracklength);

% Loop discriminators
trackResults.dllDiscr       = inf(1, tracklength);
trackResults.blksize = zeros(1, tracklength);

%C/No
trackResults.CNo.VSMValue = ...
    zeros(1,floor(tracklength/settings.CNo.VSMinterval));
trackResults.CNo.VSMIndex = ...
    zeros(1,floor(tracklength/settings.CNo.VSMinterval));

trackResults.CNo.PRMValue=0; %To avoid error message when
trackResults.CNo.PRMIndex=0; %tracking window is closed before completion.

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%% init by zsh
NumChan=length(activeChnList);%settings.numberOfChannels;
stateno=17;%number of states

%process noise var-covariance matrix
Qw=diag([diag(1e0*eye(3))',diag(1e-3*eye(3))',1*diag(1e-2*eye(3))',1*diag(1e-8*eye(3))',1*diag(1e-8*eye(3))',1e-6,1e-1]);

%measurement noise var-covariance matrix
R(1:NumChan,1:NumChan)=1500*eye(NumChan);%3000*eye(2*NumChan); % dimension=number of available channels
R(NumChan+1:2*NumChan,NumChan+1:2*NumChan)=9e2*eye(NumChan);%8e1*eye(NumChan);%

% initial estimation error var-covairance matrix
P0=diag([1e0,1e0,1e0,1e-1,1e-1,1e-1,1*diag(1e-10*eye(3))',1*diag(1e-10*eye(3))',1*diag(1e-10*eye(3))',1,1e-8]);

%initialize measurement matrix
H=zeros(2*NumChan,stateno);

%initialize measurement vector
Z=zeros(2*NumChan,1);
% states of Kalman filter initialization
X_est=zeros(stateno,tracklength);
X0=zeros(stateno,1);
ddt0=-navSolutions.ddt(StartTime*settings.navSolRate/1000);
X0(16)=-ddt0/1000;
dt=1*X0(16);

%INS measurements
load C:\Sara\School\ASEN-6080StatOD2\Project1\inputs\SPANIMUreadings
addpath 'insfunctions'
j=1;
npts = tracklengthall/10+1;%number of points in IMU measurement dataset
%find the true position for the starttime sample point
ind1=find(INSCorrTime>=navSolutions.rxTime(StartTime*settings.navSolRate/1000),1);
lat0=INSlat(ind1)*d2r;
lon0=INSlon(ind1)*d2r;
hei0=INShei(ind1);
[pos0(1,1),pos0(1,2),pos0(1,3)]=geo2cart([lat0*r2d,0,0],[lon0*r2d,0,0], hei0, 5);
pos_kf=pos0;
%find attitude for the starttime sample point
ind0=find(GPSTime>=navSolutions.rxTime(StartTime*settings.navSolRate/1000),1);

phi=INSroll(ind0)/r2d;
theta=INSpitch(ind0)/r2d ;
psi=INShead(ind0)/r2d;
%direction cosine matrix
DCMnb=eul2dcm([phi theta psi]);
%initialize the output to be saved
est_roll_KF=zeros(1,tracklength+1);
est_pitch_KF=zeros(1,tracklength+1);
est_yaw_KF=zeros(1,tracklength+1);
est_roll_KF(1) = phi;
est_pitch_KF(1) = theta;
est_yaw_KF(1) = psi;
ve=INSve(ind0);
vn=INSvn(ind0);
vu=INSvu(ind0);
%initialize intermediate variables for INS update
[tlat,tlon,thei]=cart2geo(pos0(1,1),pos0(1,2),pos0(1,3),5);
orginllh=[tlat*d2r,tlon*d2r,thei];
est_lat=zeros(1,tracklength+1);
est_lon=zeros(1,tracklength+1);
est_height=zeros(1,tracklength+1);
est_lat(1)=orginllh(1);
est_lat(2)=orginllh(1);
est_lon(1)=orginllh(2);
est_height(1)=orginllh(3);
height = orginllh(3); 

heightold = height;
veold = ve;
vnold = vn;
vuold=vu;
vel_l(1,:) = [veold vnold vu];
velenu=vel_l(1,:);
velold = [ve, vn, vu];
latold = orginllh(1);
est_DCMbn = DCMnb';
est_DCMbn_KF = est_DCMbn;
slat=sin(est_lat(1));   clat=cos(est_lat(1));
slon=sin(est_lon(1));   clon=cos(est_lon(1));
est_DCMel=[-slon -slat*clon clat*clon
        clon -slat*slon clat*slon
        0    clat       slat]';
est_DCMel_KF=est_DCMel;
omega_e=7.292115e-5;
omega_ie_E=[0 0 omega_e]';
C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU
%transform IMU raw measurements to delta velocity and delta theta in 1ms
ind0=find(INStimetag>=navSolutions.rxTime(StartTime*settings.navSolRate/1000),1);
raw_dv=0.001*[INSaccy(ind0:ind0+npts-1)-0*INSaccby(ind1:ind1+npts-1),INSaccx(ind0:ind0+npts-1)-0*INSaccbx(ind1:ind1+npts-1),-INSaccz(ind0:ind0+npts-1)+0*INSaccbz(ind1:ind1+npts-1)];
raw_dtheta=0.001*d2r*[INSgyroy(ind0:ind0+npts-1)-0*INSgyrody(ind1:ind1+npts-1),INSgyrox(ind0:ind0+npts-1)-0*INSgyrodx(ind1:ind1+npts-1),-INSgyroz(ind0:ind0+npts-1)+0*INSgyrodz(ind1:ind1+npts-1)];
%adaptive filtering window for the Kalman filter measurements
cnt=1;
lastn=500;
mat1=zeros(NumChan,lastn);
mat2=zeros(NumChan,lastn);
%initialize intermediate variables
dSv=zeros(1,NumChan);
dPlos=zeros(1,NumChan);
Vs=zeros(1,NumChan);
dVlos=zeros(1,NumChan);
carrError=zeros(1,NumChan);
carrErrorold=zeros(1,NumChan);
codeError=zeros(1,NumChan);
codeErrorold=zeros(1,NumChan);
carrFreq=zeros(1,NumChan);
codeFreq=zeros(1,NumChan);
initsample=zeros(1,NumChan);
initsampleforcode=zeros(1,NumChan);
codePhase=zeros(1,NumChan);
codePhaseStep=zeros(1,NumChan);
remCarrPhase=zeros(1,NumChan);
vsmCnt=zeros(1,NumChan);
satPosenu = zeros(3,NumChan);
satPosenu0 = zeros(3,NumChan);
satVelenu = zeros(3,NumChan);
%initialize code, frequency, transmit time
for channelNr = 1:NumChan%settings.numberOfChannels

    % Only process if PRN is non zero (acquisition was successful)

        trackResults(activeChnList(channelNr)).PRN     = trackRes(1,activeChnList(channelNr)).PRN;
        carrFreq(1,channelNr)=trackRes(1,activeChnList(channelNr)).carrFreq(StartTime);
        codeFreq(1,channelNr)=trackRes(1,activeChnList(channelNr)).codeFreq(StartTime);
        initsample(1,channelNr)=ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(StartTime));
        initsampleforcode(1,channelNr)=ceil(trackRes(1,activeChnList(channelNr)).absoluteSample(StartTime-1));
        codePhase(1,channelNr)=(initsampleforcode(1,channelNr)-trackRes(1,activeChnList(channelNr)).absoluteSample(StartTime-1))...
            /settings.samplingFreq*codeFreq(1,channelNr);
        codePhaseStep(1,channelNr) = codeFreq(1,channelNr) / settings.samplingFreq;
        tTime=...
            findTransTime(initsample(channelNr),activeChnList,svTimeTable,trackRes);
        transmitTime(activeChnList(channelNr))=tTime(activeChnList(channelNr));
        remCarrPhase(1,channelNr)=0;
 
        %C/No computation
        vsmCnt(channelNr)  = 0;
        % Get a vector with the C/A code sampled 1x/chip
        caCode0 = generateCAcode(trackRes(1,activeChnList(channelNr)).PRN);
        % Then make it possible to do early and late versions
        caCode(channelNr,:) =[caCode0(1023) caCode0 caCode0(1)];
end % for channelNr
transmitTime0=transmitTime;
blksize0=pdi*settings.samplingFreq*ones(1,NumChan);
timet=(0:2*blksize0)/settings.samplingFreq;
blksize=pdi*settings.samplingFreq*ones(1,NumChan);
samplepos=initsample;
mininit=min(initsample);
minpos=mininit;
IP1=zeros(1,NumChan);
QP1=zeros(1,NumChan);
% prepare sin and cos tables for faster carrier replication - zsh
% qstep=2*pi/10000;
% sintable=sin(0:qstep:2*pi-qstep);
% costable=cos(0:qstep:2*pi-qstep);
        
clear navSolutions
clear GPSTime INSCorrTime INSaccbx INSaccby INSaccbz INSaccx INSaccy INSaccz
clear INSarx INSary INSarz INSgyrodx INSgyrody INSgyrodz
clear INSgyrox INSgyroy INSgyroz INShead INSpitch INSroll INStimetag INSve INSvu INSvn
%% Start processing channels ==============================================
for sectcnt=1:sectno
    i=2;
        %=== Process the number of specified code periods =================
        for loopCnt =  1:tracklength
            %calculate transmit time for every millisecond
            transmitTime(activeChnList)=transmitTime(activeChnList)-...
                (-dt)/settings.c+blksize/settings.samplingFreq-(dSv-dPlos)/settings.c;
            if rem(loopCnt-1,100)==0%update sv position per 100ms to reduce calculation
                [satPositionsall, satClkCorrall] = satpos1([transmitTime(transmitTime>0),transmitTime0(transmitTime0>0)], ...
                        [trackRes(activeChnList).PRN,trackRes(activeChnList).PRN],eph);
                satPositions=satPositionsall(:,1:NumChan);
                satPositions0=satPositionsall(:,end/2+1:end);
            else
                satPositions(1:3,:)=satPositions(1:3,:)+kmt*satPositions(4:6,:);
                satPositions0(1:3,:)=satPositions0(1:3,:)+kmt*satPositions0(4:6,:);
            end
            %read dataset segment into memory
            fseek(fid, ...
                dataAdaptCoeff*(0*settings.skipNumberOfSamples + ...
                minpos-1), ...
                'bof');
            [rawSignal0, samplesRead0] = fread(fid, ...
                    dataAdaptCoeff*(max(samplepos)-minpos+max(blksize)), settings.dataType);

            for m=1:NumChan
                rawSignal= rawSignal0((samplepos(m)-minpos)*dataAdaptCoeff+1:(samplepos(m)-minpos+blksize(1,m))*dataAdaptCoeff)';

                if (dataAdaptCoeff==2)
                    rawSignal1=rawSignal(1:2:end);
                    rawSignal2=rawSignal(2:2:end);
                    rawSignal = rawSignal1 + 1i.* rawSignal2;  %transpose vector
                end

                if rem(loopCnt-1,100)==0%transform sv pos to enu per 100ms to reduce calculation
                    satPosenu(1:3,m)=est_DCMel_KF*(satPositions(1:3,m)-pos_kf');
                    satPosenu0(1:3,m)=est_DCMel_KF*(satPositions0(1:3,m)-pos_kf');
                    satVelenu(1:3,m)=(satPosenu(1:3,m)-satPosenu0(1:3,m))/kmt;
                else
                    dvtmp=(satVelenu(1:3,m)+vel_l(i-1,1:3)')*kmt;
                    satPosenu(1:3,m)=satPosenu(1:3,m)+dvtmp;
                    satPosenu0(1:3,m)=satPosenu0(1:3,m)+dvtmp;
                end
                %calculate LOS
                le=satPosenu(1,m);
                ln=satPosenu(2,m);
                lu=satPosenu(3,m);
                norm_a=sqrt(le*le+ln*ln+lu*lu);
                a=[le;ln;lu]/norm_a;
                %form measurement matrix
                H(m,:)=[-a(1),-a(2),-a(3),zeros(1,12),-1,0];
                H(m+NumChan,:)=[zeros(1,3),+a(1),+a(2),+a(3),zeros(1,9),0,1];
                
                dSv(m)=(satPosenu(1:3,m)-satPosenu0(1:3,m))'*a;%sv displacement projection on LOS
                Vs(m)=(satVelenu(1:3,m)'-velenu)*a;%relative velocity between sv and user on LOS
               
                dPlos(m)=(kmt*velenu)*a;%the user displacement between current and previous epoch
                dVlos(m)=X0(4:6)'*a;%estimated user velocity error on LOS
                %update code frequency and phase
                codeFreq(1,m)=settings.codeFreqBasis*(1-(ddt0+Vs(m))/settings.c);%ddt0/settings.c-Vs(m)/settings.c);
                codePhaseStep(1,m) = codeFreq(1,m) / settings.samplingFreq;                
                %correct previous codephase with estimated position error
                codePhase(1,m) =codePhase(1,m) + (dt+X0(1:3)'*a)/settings.c*codeFreq(1,m);%dt/settings.c*codeFreq(1,m)+X0(1:3)'*a/settings.c*codeFreq(1,m);
                %generate current 1ms code phase
                codePhase(1,m) =codePhase(1,m) -...
                    (dSv(m)-dPlos(m))/settings.c*codeFreq(1,m)+(blksize(1,m)-blksize0(1,m)).*codePhaseStep(1,m);
                %generate early code
%                 tcode       = (codePhase(1,m)-earlyLateSpc) : ...
%                     codePhaseStep(1,m) : ...
%                     ((blksize(1,m)-1)*codePhaseStep(1,m)+codePhase(1,m)-earlyLateSpc);
%                 tcode2      = ceil(tcode+ 1);
%                 earlyCode   = caCode(m,tcode2);

                % Define index into late code vector
%                 tcode       = (codePhase(1,m)+earlyLateSpc) : ...
%                     codePhaseStep(1,m) : ...
%                     ((blksize(1,m)-1)*codePhaseStep(1,m)+codePhase(1,m)+earlyLateSpc);
%                 tcode2      = ceil(tcode+ 1) ;
%                 lateCode    = caCode(m,tcode2);

                % Define index into prompt code vector
                tcode       = codePhase(1,m) : ...
                    codePhaseStep(1,m) : ...
                    ((blksize(1,m)-1)*codePhaseStep(1,m)+codePhase(1,m));
                tcode2      = ceil(tcode+ 1) ;
                promptCode  = caCode(m,tcode2);
                ELspcpnt=ceil(earlyLateSpc/codePhaseStep(1,m));
                %about 4% loss compared to the commented lines, but faster - zsh
                earlyCode=[promptCode(end-ELspcpnt+2:end),promptCode(1:end-ELspcpnt+1)];
                lateCode=[promptCode(ELspcpnt:end),promptCode(1:ELspcpnt-1)];

                %% Generate the carrier frequency to mix the signal to baseband -----------
                carrFreq(1,m)=settings.IF+(ddt0+Vs(m))/settings.c*1575.42e6;
                time=timet(1:blksize(1,m)+1);

                % Get the argument to sin/cos functions
                trigarg = ((carrFreq(1,m) * 2.0 * pi) .* time) + remCarrPhase(1,m);
                remCarrPhase(1,m) = mod(trigarg(blksize(1,m)+1), (2 * pi));
%                 indremCP=floor(mod(trigarg(1:blksize(1,m)),(2*pi))/qstep)+1;
                %using table to generate carrier replicas
%                 carrsigsin=sintable(indremCP);
%                 carrsigcos=costable(indremCP);

                % Finally compute the signal to mix the collected data to bandband
                carrsig = exp(1i .* trigarg(1:blksize(1,m)));%so time-consuming, to be optimized...
%                 carrsig=carrsigcos+1i*carrsigsin;

                %% Generate the six standard accumulated values ---------------------------
                % First mix to baseband
                qBasebandSignal = real(carrsig .* rawSignal);
                iBasebandSignal = imag(carrsig .* rawSignal);

                % Now get early, late, and prompt values for each
                I_E = earlyCode  * iBasebandSignal';%more time-efficient than previous sum function
                Q_E = earlyCode  * qBasebandSignal';
                I_P = promptCode * iBasebandSignal';
                Q_P = promptCode * qBasebandSignal';
                I_L = lateCode   * iBasebandSignal';
                Q_L = lateCode   * qBasebandSignal';

                
                %% Find PLL error and update carrier NCO ----------------------------------

                % Implement carrier loop discriminator (phase detector)

                if (loopCnt==1 && sectcnt==1)
                    IP1(1,m)=I_P;
                    QP1(1,m)=Q_P;
                    carrErrorold(1,m)=carrError(1,m);
                    carrError(1,m)=0;
                else
                    dot=IP1(1,m)*I_P+QP1(1,m)*Q_P;
                    cross=IP1(1,m)*Q_P-I_P*QP1(1,m);
                    % frequency discriminator
                    carrErrorold(1,m)=carrError(1,m);
                    carrError(1,m) = cross*sign(dot)/(2*pi*(I_P*I_P+Q_P*Q_P));
                    IP1(1,m)=I_P;
                    QP1(1,m)=Q_P;
                end
                
                % Implement carrier loop filter and generate NCO command
                Z(m+NumChan,1)=(carrErrorold(1,m)+(carrError(1,m)-carrErrorold(1,m))...
                    /blksize(1,m)*(blksize(1,m)-rem((samplepos(1,m)-minpos),blksize(1,m))))/pdi/1575.42e6*settings.c;

                trackResults(activeChnList(m)).carrFreq(loopCnt) = carrFreq(1,m);

                %% Find DLL error and update code NCO -------------------------------------
                codeErrorold(1,m)=codeError(1,m);
                Et=(I_E * I_E + Q_E * Q_E);
                Lt=(I_L * I_L + Q_L * Q_L);
%                 codeError(1,m) = ((I_E * I_E + Q_E * Q_E) - (I_L * I_L + Q_L * Q_L)) / ...
%                     2/((I_E * I_E + Q_E * Q_E) + (I_L * I_L + Q_L * Q_L));
                codeError(1,m) = (Et - Lt)/2/(Et + Lt);

                Z(m,1)=(codeErrorold(1,m)+(codeError(1,m)-codeErrorold(1,m))...
                    /blksize(1,m)*(blksize(1,m)-rem((samplepos(1,m)-minpos),blksize(1,m))))/codeFreq(1,m)*settings.c;
                %% Record various measures to show in postprocessing ----------------------
%comment the following lines for speed - zsh
%                 trackResults(activeChnList(m)).codeFreq(loopCnt) = codeFreq(1,m);
%                 trackResults(activeChnList(m)).remCodePhase(loopCnt)=codePhase(1,m);
%                 trackResults(activeChnList(m)).dllDiscr(loopCnt)       = codeError(1,m);
%                 trackResults(activeChnList(m)).I_E(loopCnt) = I_E;
%                 trackResults(activeChnList(m)).I_P(loopCnt) = I_P;
%                 trackResults(activeChnList(m)).I_L(loopCnt) = I_L;
%                 trackResults(activeChnList(m)).Q_E(loopCnt) = Q_E;
%                 trackResults(activeChnList(m)).Q_P(loopCnt) = Q_P;
%                 trackResults(activeChnList(m)).Q_L(loopCnt) = Q_L;

%                 if (settings.CNo.enableVSM==1)
%                     if (rem(loopCnt,settings.CNo.VSMinterval)==0)
%                         vsmCnt(m)=vsmCnt(m)+1;
%                         CNoValue=CNoVSM(trackResults(activeChnList(m)).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
%                             trackResults(activeChnList(m)).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
%                         trackResults(activeChnList(m)).CNo.VSMValue(vsmCnt(m))=CNoValue;
%                         trackResults(activeChnList(m)).CNo.VSMIndex(vsmCnt(m))=loopCnt;
%                     end
%                 end
%                 trackResults(activeChnList(m)).blksize(loopCnt)  = blksize(1,m);
            end

            if (rem(loopCnt-1,10)==0)%to sync 1ms loop and 10ms INS
                j=j+1;
            end
            %% INS strapdown update
            tupd = kmt;%update interval
            %update dcm from epoch k to epoch k+1 using gyro delta angles
            est_DCMbn = est_DCMbn_KF*calDCM(raw_dtheta(j,1:3)-1*X0(13:15)'*tupd);
            %rotation rate of enu relative to the ecef, expressed in enu (rad/s)
            omega_el_L = llangrate(latold,veold,vnold,heightold);
            omega_ie_L = est_DCMel_KF*omega_ie_E;% earth rotation rate relative to inertial frame expressed in ENU
            omega_il_L = omega_ie_L + omega_el_L;% enu rotation relative to inertial, expressed in ENU
            DCM_ll_I = calDCM(-omega_il_L*tupd);%enu dcm between epochs
            est_DCMbn = C*(DCM_ll_I*(C*est_DCMbn)); %estimated DCM b to n, taking local-level frame rotation into account
            est_delv_b = raw_dv(j,1:3)-1*X0(10:12)'*tupd;       % extract delta-V for current point in time
            del_Vl = C*(est_DCMbn*est_delv_b');
%             if (rem(loopCnt-1,10)==0)%reduce computation
                % gravity changes slower, can be computed less
                g0 = gravity(latold,heightold);
%             end
            
            est_DCMel=calDCM(-omega_el_L*tupd)*est_DCMel_KF;
            omega_ie_L = est_DCMel*omega_ie_E;
            as=antisymm([veold,vnold,vuold]);
            vtmp=velold'+del_Vl+(as*(omega_el_L+2*omega_ie_L)+[0, 0, -g0]')*tupd;
            vel_l(i,:) = vtmp';

            est_height(i) = est_height(i-1)+vel_l(i,3)*tupd;
            heightold = est_height(i);
            radiusa=localradius(est_lat(i));
            accel_L = del_Vl/tupd;
            
            est_lat(i) = asin(est_DCMel(3,3));
            est_lon(i) = atan2(est_DCMel(3,2),est_DCMel(3,1));
            slat=sin(est_lat(i));
            clat=cos(est_lat(i));
            ve=vel_l(i,1);
            vn=vel_l(i,2);
            vu=vel_l(i,3);
            rh=radiusa+est_height(i);
            
            %error propagation coefficient matrix
            F11 = zeros(3,3);    % L-frame
            F11(1,2)=ve*slat/rh/clat/clat;
            F11(1,3)=ve/rh/rh/clat;%N and E
            F11(2,3)=vn/rh/rh;%N and U
            
            F12 = eye(3);
            F13 = zeros(3,3);
            
            F21=zeros(3,3);
            F21(1,2)=-ve*(2*omega_e*clat+ve/rh/clat/clat);
            F21(1,3)=(ve*ve*slat/clat+vn*vu)/rh/rh;
            F21(2,2)=2*omega_e*(vn*clat+vu*slat)+vn*ve/rh/clat/clat;
            F21(2,3)=-ve/rh/rh*(vn*slat/clat-vu);
            F21(3,2)=-2*omega_e*ve*slat;
            F21(3,3)=(ve*ve+vn*vn)/rh/rh;
            
            F22 = zeros(3,3);
            F22(1,1)=(vn*slat/clat-vu)/rh;
            F22(1,2)=2*omega_e*slat+ve*slat/clat/rh;
            F22(1,3)=-(2*omega_e*clat+ve/rh);
            F22(2,1)=-2*(omega_e*slat+ve*slat/clat/rh);
            F22(2,2)=-vu/rh;
            F22(2,3)=-vn/rh;
            F22(3,1)=2*(omega_e*clat+ve/rh);
            F22(3,2)=2*vn/rh;
            
            F23 = antisymm(accel_L);
            F31 = zeros(3,3);
            F32 = zeros(3,3);
            F33 = (-1)*antisymm(omega_il_L);

            F = zeros(17,17);
            F(1:9,1:9) = [F11 F12 F13; F21 F22 F23; F31 F32 F33]; 
            F(10:12,10:12) = zeros(3,3);
            F(4:6,10:12) = C*est_DCMbn;
            F(13:15,13:15) = zeros(3,3);
            F(7:9,13:15) = C*est_DCMbn;
            F(16,17) = 1;
            %discretize F matrix, assuming all states are irrelavent to
            %themselves at last epoch, except for clock bias and drift
            F=F*kmt+eye(stateno);
            F(1,1)=0;
            F(2,2)=0;
            F(3,3)=0;
            F(4,4)=0;
            F(5,5)=0;
            F(6,6)=0;
            F(7,7)=0;
            F(8,8)=0;
            F(9,9)=0;
            F(16,16)=1;
            F(17,17)=0;
            %Kalman filter
            P=F*P0*F'+Qw;
            K=P*H'/(H*P*H'+R);
            P0=(eye(stateno)-K*H)*P;
            X_next=F*X0;

            beta=(Z(:,1)-H*X_next(:,1));
            alpha=K*beta;

            X_est(:,loopCnt)=X_next(:,1)+alpha;%
            X0=X_est(:,loopCnt);

         %% adaptive filter
            res=beta;
            mat1(:,cnt)=res(1:NumChan);%code
            mat2(:,cnt)=res(NumChan+1:end);%carrier

            cnt=cnt+1;
            if cnt==lastn
                cnt=1;
            end
            if (loopCnt>lastn ||sectcnt>1)
                Cres1=mat1*mat1'/lastn;
                Cres2=mat2*mat2'/lastn;
                tmpR=diag([diag(Cres1);diag(Cres2)]+diag(H*P0*H'));
                R(1:NumChan,1:NumChan)=tmpR(1:NumChan,1:NumChan);
                R(NumChan+1:2*NumChan,NumChan+1:2*NumChan)=tmpR(NumChan+1:2*NumChan,NumChan+1:2*NumChan);
            end
            %correct errors
            est_lat(i)=est_lat(i)+X0(2)/radiusa;
            est_lon(i)=est_lon(i)+X0(1)/radiusa/cos(est_lat(i));
            est_height(i)=est_height(i)+1*X0(3);
            theta(1,1) = -X_est(2,loopCnt)/radiusa;
            theta(2,1) = X_est(1,loopCnt)/radiusa;
            theta(3,1) = tan(est_lat(i))*theta(2);

            psi=X0(7:9);
            phi_angle=psi+theta;%total attitude error
            slat=sin(est_lat(i));
            clat=cos(est_lat(i));
            slon=sin(est_lon(i));
            clon=cos(est_lon(i));
            est_DCMel_KF=[-slon -slat*clon clat*clon
                    clon -slat*slon clat*slon
                    0    clat       slat]';
            latold = est_lat(i);
            est_DCMbn_KF = C*(eye(3)+antisymm(phi_angle))*C*est_DCMbn; 
            eulangle = dcm2eul(est_DCMbn_KF);
            est_roll_KF(i) = eulangle(1);
            est_pitch_KF(i) = eulangle(2);
            est_yaw_KF(i) = eulangle(3);

            vel_l(i,:)=vel_l(i,:)-X0(4:6)';
            veold=vel_l(i,1);
            vnold=vel_l(i,2);
            vuold=vel_l(i,3);
            velold = vel_l(i,:);
            dt=X0(16);
            ddt=X0(17);

            ddt0=ddt0+ddt;
            [pos_kf(1,1),pos_kf(1,2),pos_kf(1,3)]=geo2cart([est_lat(i)*r2d,0,0],[est_lon(i)*r2d,0,0], est_height(i), 5);
            vel_kf=(est_DCMel_KF'*vel_l(i,:)')';
            velenu=vel_l(i,:);
            
            pos_kf=vel_kf*kmt+pos_kf;%estimate next pos in ecef
            
            transmitTime0=transmitTime;
            blksize = ceil((settings.codeLength-codePhase) ./ codePhaseStep);

            samplepos=samplepos+blksize;
            minpos=min(samplepos);
            i=i+1;
        end % for loopCnt
%         disp('   Saving VLL results to file "deepresult.mat"')
        fname=['deepresult',num2str(StartTime/1000+0*tracklength/1000*(sectcnt-1)),...
            '~',num2str(StartTime/1000+tracklength/1000)];
        save(fname, ...
             'settings', 'StartTime','trackResults', 'tracklength','activeChnList',...
            'P0','samplepos','transmitTime','X_est',...
            'est_lon','est_lat','est_height','est_roll_KF','est_pitch_KF','est_yaw_KF','ind1','Qw','activeChnList','npts','F');
        %prepare for next tracking section
        StartTime=StartTime+tracklength;
        vel_l(1,:)=vel_l(i-1,:);
        est_height(1)=est_height(i-1);
        fprintf('time consumed:%f, finish time:%s\n',toc,datestr(now));
end
% disp('tracking finished')
fprintf('tracking completed. time consumed:%f, finish time:%s\n',toc,datestr(now));
plotDeep(begintime,endtime);
toc