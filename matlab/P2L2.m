function [Mxyev, Mxyiv, vb, Mzev, Mziv] = P2L2(vars,fdv)
%Format: function [Mxyev, Mxyiv, vb, Mev, Miv] = P2LBv4(vars,fdv)
%
%input:
    % vars:                 values for parms have been or will be fit
    % fdv = "forward variables" containing info about the measurement
    %    fdv.fitvars        Names of parms that need to be fit
    %    fdv.knowns         Names of parms that are known
    %    fdv.knownvals      Values of known parameters
    %    fdv.FlipAngles     Excitation angle list [2 x NFlips] in degrees
    %    fdv.TR:            Repetition time [2 x NFlips]
    %    fdv.taxis:         time axis for the data [1 x NFlips]
    %    fdv.data           [2 x NTP] observed [pyr;lac] data
    %    fdv.UseVIF         =1 if VIFs have been measured
    %    fdv.VIFP           **Mz** VIF for pyruvate [1 x NFlips]
    %    fdv.VIFL           **Mz** VIF for lactate [1 x NFlips] (usually zeros)
    %    fdv.verbose        =1 to plot results of this script
%    
%output:
    % Mxyev                 [2 x NTP] extravascular observed [Pyr;Lac] 
    % Mxyiv                 [2 x NTP] intravascular observed [Pyr;Lac]
    % vb                    vascular blood vol fraction
    % Mzev                   [2 x NTP] extravascular longitudinal [Pyr;Lac]
    % Mziv                   [2 x NTP] intravascular observed [Pyr;Lac]
%
% Note Mzev,Mziv reflects z-directred magnetization before the excitation pulse ,  
%  while Mxy reflects transverse magnetization due to pulse.
%
% Note total observed signal is [Pyr;Lac] = (1-vb)*Mxyev + vb*Mxyiv
%

    eps=1e-6;
       
    %Unpack variables:
    for ii=1:length(vars)
        x.(fdv.fitvars{ii})=vars(ii);
    end
    for ii=1:length(fdv.knowns)
        x.(fdv.knowns{ii})=fdv.knownvals(ii);
    end
    vb=x.vb;
    kvedve=x.kve/(1-vb);

    %Assume for now that the VIF is given:
    if fdv.UseVIF
        ff=[x.VIFScale*fdv.VIFP;
            x.VIFScale*fdv.VIFL];
    else 
        ff=[x.VIFScale*gampdf(fdv.taxis,x.Gam1,x.Gam2);
            zeros(1,fdv.NFlips)];
    end
    %IV magnetization at each segment:
    MivSeg =ff;
    MxyivSeg=ff.*sind(fdv.FlipAngle);
    
    %Initialize return variables:
    Mxyev=zeros(2,fdv.ntp);
    Mzev =zeros(2,fdv.ntp);
    Mxyiv=zeros(2,fdv.ntp);
    Mziv =zeros(2,fdv.ntp);

    %Manage trivial case of vb=0: no signal, or vb=1: all vessel
    if (vb<eps) 
        %All return variables have been initialized to zero, above.
        return
    elseif (vb>(1-eps))
        % If segmented acquisition with NFlips>NTP, where VIF has been
        % measured or interpolated to each excitation
        % Combine segments: ** NOTE HERE ASSUMES EQUAL WEIGHTING **
        for ii=1:fdv.ntp 
            Mxyiv(:,ii)=mean(MxyivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
            Mziv(:,ii)=mean(MivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
        end
        return
    end
    
    
    %Calculate nontrivial result
    % [Pev;Lev] = [Pev0;Lev0]*exp(At) +kpl*integral(0,t:exp(A*(t-tau))*[Piv(tau);Liv(tau)]*dtau)
    % A is 2x2
    a11=-((kvedve) + x.kpl + (1/x.T1Pyr));
    a12=x.klp;
    a21=x.kpl;
    a22=-((kvedve) + x.klp + (1/x.T1Lac));
    A = [a11 a12 ; a21 a22];
    % Diagonlize to permit matrix integral.  
    [P,D]=eig(A);
    % D is diagonlized matrix (2x2)
    % P is matrix of eigenvectors 
    %   --> A*P=P*D
    %   --> A=P*D\P
    dD=diag(D); % get diagonal values (1x2)
    
    %Calculate signal evolution each TR:
    
    MzevSegIC=[x.P0; x.L0]; %Set Initial Conditions
    for ii=1:fdv.ntp
        MxyevSeg=zeros(2,fdv.NSeg);
        MzevSeg =zeros(2,fdv.NSeg);
        for jj=1:fdv.NSeg
            iiseg=(ii-1)*(fdv.NSeg)+jj; % Works when NSeg=1 too
            TR=fdv.TR(iiseg);
            %
            %First account for EV signal already present and its evolution
            %
            %Longitudinal magnetization available at the start of each segment:
            MzevSeg(:,jj)=MzevSegIC;
            %Signal observed at each excitation, at start of segment:
            MxyevSeg(:,jj)=MzevSegIC.*sind(fdv.FlipAngle(:,iiseg));
            %At the end of this TR, Mz evolves to contribute to IC for next:
            if iiseg<fdv.NFlips % Don't need to calc after last datapoint
                MzevSeg1=(exp(dD*TR)).*(P\(MzevSegIC.*cosd(fdv.FlipAngle(:,iiseg))));
                %
                %Now calculate new spins flowing into the system
                %
                %Assume piecewise linear VIF. Diagonalize:
                dff1=P\ff(:,iiseg);    % Diag'd VIF @ start of TR
                dff2=P\ff(:,iiseg+1);  % Diag'd VIF @ end of TR
                %Get slope and y-intercept for diagonolized forcing function:
                b=dff1;
                m=(dff2-dff1)/TR;
                %At the end of this TR, inflowing spins will lead to:
                %MevSeg2a=exp(dD*TR).*(-b./dD).*(exp(-dD*TR)-1);
                %MevSeg2b=exp(dD*TR).*m.*(((-TR./dD)-(1./dD./dD)).*exp(-dD*TR)+(1./dD./dD));
                %slightly more efficient way to calculate:
                MevSeg2a=(-b./dD).*(1-exp(dD*TR));
                MevSeg2b=m.*(((-TR./dD)-(1./dD./dD))+exp(dD*TR).*(1./dD./dD));
                %Total signal at end of TR equals IC for next TR:
                MzevSegIC = P*(MzevSeg1 + kvedve*(MevSeg2a+MevSeg2b));
               
            end
        end
        if (fdv.NSeg>1) %manage signal from segmented acquisitions
            %sizemxyev=size(mean(MxyevSeg))
            Mxyev(:,ii)=mean(MxyevSeg,2);
            Mzev(:,ii) =mean(MzevSeg,2);
            Mxyiv(:,ii)=mean(MxyivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
            Mziv(:,ii) =mean(MivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
        else
            Mxyev(:,ii)=MxyevSeg;
            Mzev(:,ii) =MzevSeg;
            Mxyiv(:,ii)=MxyivSeg(:,ii);
            Mziv(:,ii) =MivSeg(:,ii);
        end            
    end

if fdv.verbose
    figure(99)
    plot(fdv.taxis(1:fdv.NSeg:end),Mxyev(1,:),'g-',...
        fdv.taxis(1:fdv.NSeg:end),Mxyev(2,:),'b-')
    legend('PyrEV','LacEV')
end    