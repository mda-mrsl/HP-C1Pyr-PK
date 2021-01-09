%% Set up model system:
clear

%Unknowns:
fdv.fitvars={'kpl' 'VIFScale'};

%Nuisance parameters; known or estimated eslewhere.
%Noise terms estimate mean value of (magnitude) noise with nonzero mean;
% set noise terms to -Inf when real or complex data is being used.
fdv.knowns={'kve','vb','T1Pyr','T1Lac','klp','Gam1','Gam2','tdel','P0','L0'};
fdv.knownvals=[0.02, 0.05, 43, 33, 0, 2.8, 4.5, 10, 0, 0];



%Describe acquisition scheme
fdv.ntp=64; % Number of timepoints
fdv.NSeg=1; % Segments per timepoint
fdv.NFlips=(fdv.ntp)*(fdv.NSeg); % Total number of excitations

%Describe temporal sampling scheme
TR=2; % constant repetition time
fdv.TR=ones(1,fdv.NFlips)*TR;
fdv.taxis=cumsum(fdv.TR)-fdv.TR(1);

%Describe excitation scheme
fdv.FlipAngle=20*ones(2,fdv.NFlips);

%Describe vascular inmput function
fdv.UseVIF=1;
fdv.VIFP=gampdf(fdv.taxis-6,2.8,4.5);
fdv.VIFL=zeros(1,fdv.NFlips);

%Placeholder for data to be fit
fdv.data=ones(2,fdv.ntp);
fdv.Name='HP Model Demo';

%Other miscellaneous
fdv.verbose=0;

%% Generate synthetic data:

parms=[0.1 1000]; %values for variables in fdv.fitvars
[EV, IV, vb]=P2L2(parms,fdv);

fdv.data=IV*vb+EV*(1-vb);

%resid=P2LBv4Err(parms,fdv);

figure(1);clf
plot(fdv.taxis,fdv.data(1,:),'g-',...
    fdv.taxis,fdv.data(2,:),'b-');
legend('Total Pyr','Total Lac')

%% Fit synthetic data and try to recover parms defined above

nruns=100;
jbopts=optimset('display','off');
LB   =[0 0];
UB   =[1 Inf];

bestresid=Inf;
bestfits = zeros(1,length(parms));
for ii=1:nruns
    %ii
    Guess=(parms*10).*rand(1,2);
    [fits,resid] = lsqnonlin(@(x) P2L2Err(x,fdv),Guess,LB,UB,jbopts);
    if resid<bestresid
        bestresid=resid;
        bestfits=fits;
    end
end

fprintf('Kpla in = %5.3f; Kpla out = %5.3f\n',parms(1),bestfits(1));

%% Plot

fh=P2L2Plot(bestfits,fdv);


