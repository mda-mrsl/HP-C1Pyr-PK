function [Mxy, Mz] = P2L1(vars,fdv)
%Format: function [Mxyev, Mxyiv, vb, Mev, Miv] = P2LBv4(vars,fdv)
%
%input:
    % vars:                 values for parms have been or will be fit
    % fdv = "forward variables" containing info about the measurement
    %    fdv.fitvars        Names of parms that need to be fit
    %    fdv.knowns         Names of parms that are known
    %    fdv.knownvals      Values of known parameters
    %    fdv.FlipAngle      Excitation angle list [2 x NFlips] in degrees
    %    fdv.TR:            Repetition time [2 x NFlips]
    %    fdv.taxis:         time axis for the data [1 x NFlips]
    %    fdv.data           [2 x NTP] observed [pyr;lac] data
    %    fdv.VIFP           *** NOT USED FOR THIS MODEL ***
    %    fdv.verbose        =1 to plot results of this script
 
%output:
    % Mxy                   [2 x NTP] Mxy observed [Pyr;Lac] 
    % Mz                    [2 x NTP] Mz [Pyr;Lac]
%
% Note Mz reflects magnetization before excitation pulse,  
%  while Mxy reflects transverse magnetization due to pulse.
%
% Note total observed signal is [Pyr;Lac]
%

    eps=1e-6;
       
    %Unpack variables:
    for ii=1:length(vars)
        x.(fdv.fitvars{ii})=vars(ii);
    end
    for ii=1:length(fdv.knowns)
        x.(fdv.knowns{ii})=fdv.knownvals(ii);
    end

    %Data from each excitation:
    PxySeg  = fdv.data(1,:);
    PzSeg   = PxySeg./sind(fdv.FlipAngle(1,:));
    
    %Initialize return variables:
    Mxy  = zeros(2,fdv.ntp);
    Mz   = zeros(2,fdv.ntp);

    %Set up calculation
    % Lac = Lac0*exp(At) + kpl*integral(0,t:exp(A*(t-tau))*[Pyr(tau)]*dtau)
    % Here, A is a single value:
    A = -(x.klp + (1/x.T1Lac));
    
    %Calculate signal evolution each TR:
    
    LzSegIC = x.L0; %Set Initial Conditions
    for ii=1:fdv.ntp
        LxySeg=zeros(1,fdv.NSeg);
        LzSeg =zeros(1,fdv.NSeg);
        for jj=1:fdv.NSeg
            iiseg=(ii-1)*(fdv.NSeg)+jj; % Works when NSeg=1 too
            TR=fdv.TR(iiseg);
            %
            %First account for EV signal already present and its evolution
            %
            %Longitudinal magnetization available at the start of each segment:
            LzSeg(jj)=LzSegIC;
            %Signal observed at each excitation, at start of segment:
            LxySeg(jj)=LzSegIC.*sind(fdv.FlipAngle(2,iiseg));
            %At the end of this TR, Mz evolves to contribute to IC for next:
            if iiseg<fdv.NFlips % Don't need to calc after last datapoint
                LzSeg1=exp(A*TR)*LzSegIC*cosd(fdv.FlipAngle(2,iiseg));
                %MevSeg1=(exp(dD*TR)).*(P\(MevSegIC.*cosd(fdv.FlipAngle(:,iiseg))));
                %
                %Now calculate new spins flowing into the system during this TR:
                %
                %Get slope and y-intercept for pyruvate forcing function:
                b=PzSeg(1,iiseg);
                m=(PzSeg(1,iiseg+1)-PzSeg(1,iiseg))/TR;
                %At the end of this TR, inflowing spins will lead to:
                LzSeg2= (((m/A+b)*exp(A*TR)) - (m*(TR+1/A)) - b)/A;
                %Total signal at end of TR equals IC for next TR:
                LzSegIC = LzSeg1 + (x.kpl)*LzSeg2;
            end
        end
        if (fdv.NSeg>1)
            %Note that this calculates the average signal from a segmented
            %acquisition.  May need to modify or add weighting factors to 
            %compensate for k-space coverage.
            %
            %Pyr
            Mxy(1,ii)=mean(PxySeg(1,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
            Mz(1,ii)=mean(PzSeg(1,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
            %Lac
            Mxy(2,ii)=mean(LxySeg);
            Mz(2,ii)=mean(LzSeg);
        else
            Mxy(:,ii)=[PxySeg(1,ii); LxySeg];
            Mz(:,ii)=[PzSeg(1,ii) ; LzSeg];
        end            
    end

if fdv.verbose
    figure(99)
    plot(fdv.taxis(1:fdv.NSeg:end),Mxy(1,:),'g-',...
        fdv.taxis(1:fdv.NSeg:end),Mxy(2,:),'b-')
    legend('Pyr','Lac')
    xlabel('Time (s)')
    grid on
end    