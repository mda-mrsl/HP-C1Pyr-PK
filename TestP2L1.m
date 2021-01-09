%% Generate synthetic data:
%Confirm that this algorithm works as expected:
clear

Seq=2;
if Seq==1
    % Segmented EPSI similar to UCSF:
    NTP=60;
    NSeg=8;
    TR=0.25;
    FA=acosd(cosd(20)^(1/NSeg));
elseif Seq==2
    % Single shot, typical MDA:
    NTP=60;
    NSeg=1;
    TR=2;
    FA=20;
end

fdva.fitvars={'kpl'};
fdva.knowns={'T1Lac','klp','L0'};
fdva.knownvals=[33, 0, 0];
fdva.ntp=NTP;
fdva.NSeg=NSeg;
fdva.NFlips=(fdva.ntp)*(fdva.NSeg);
fdva.TR=ones(1,fdva.NFlips)*TR;
fdva.taxis=cumsum(fdva.TR)-fdva.TR(1);
fdva.FlipAngle=FA*ones(2,fdva.NFlips);
%fdv.VIFP=gampdf(fdv.taxis-10,2.8,4.5);
fdva.data=zeros(2,fdva.NFlips);
fdva.data(1,:)=gampdf(fdva.taxis-10,2.8,4.5).*sind(fdva.FlipAngle(1,:));
fdva.Name='Test P2LAv4';
fdva.verbose=0;

% fdv=struct('fitvars',{'kpl' 'kve' 'vb' 'VIFScale'},...
%     'knowns',{'T1Pyr','T1Lac','klp','Gam1','Gam2','tdel'},...
%     'taxis',(0:2:190),'FlipAngles',20*ones(2,96),'EffBurstAngles',20*ones(2,96));

parms=[0.05]; % this is kpl only
[Mxy,Mz]=P2L1(parms,fdva);
fdva.data=Mxy;

%resid=P2LAv4Err(parms,fdv);

%% Fit

jbopts=optimset('display','off');
[fits,resid] = lsqnonlin(@(x) P2L1Err(x,fdva),0.01,0,Inf,jbopts);
fprintf('Kpla in = %5.3f; Kpla out = %5.3f\n',parms,fits);

%% Plot

fh=P2L1Plot(fits,fdva);

