function [Mxyev, Mxyiv, vols, Mzev, Mziv] = P2L3(vars,fdv)
%Format: function [Mxy, vols, Mev, Miv] = P2LBv4(vars,fdv)
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
    % Mxyev                 [4 x NTP] extravascular observed [PyrEE;LacEE;PyrIC;LacIC] 
    % Mxyiv                 [4 x NTP] IV observed [PyrIV; LacIV; 0; 0]
    % vols                  vascular blood vol fractions [vb ve vc]            
    % Mzev                  [4 x NTP] extravascular longitudinal [PyrEE;LacEE;PyrIC;LacIC]
    % Mziv                  [4 x NTP] IV longitudinal [PyrIV; LacIV; 0; 0];

%
% Note Mz reflects magnetization before excitation pulse 
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
    if ~isfield(x,{'kecl'}) %set equal if not specified
        x.kecl=x.kecp;
    end
    vb=x.vb;
    ve=(1-vb)*x.vef;
    vc=1-vb-ve;
    vols=[vb ve vc];
    R1Pyr=1/x.T1Pyr;
    R1Lac=1/x.T1Lac;
    kvedve=x.kve/ve;
    kecpdve=x.kecp/ve;
    kecpdvc=x.kecp/vc;
    kecldve=x.kecl/ve;
    kecldvc=x.kecl/vc;
    
    %Assume for now that the VIF is given:
    if fdv.UseVIF
        ff=[x.VIFScale*fdv.VIFP;
            x.VIFScale*fdv.VIFL;
            zeros(1,fdv.NFlips);
            zeros(1,fdv.NFlips)];
    else 
        ff=[x.VIFScale*gampdf(fdv.taxis,x.Gam1,x.Gam2);
            zeros(1,fdv.NFlips);
            zeros(1,fdv.NFlips);
            zeros(1,fdv.NFlips)];
    end
    
    %IV magnetization at each segment:
    MzivSeg =ff;
    MxyivSeg=ff.*sind(fdv.FlipAngle);
    
    %Initialize return variables:
    Mxyev=zeros(4,fdv.ntp);
    Mxyiv=zeros(4,fdv.ntp);
    Mzev=zeros(4,fdv.ntp);
    Mziv=zeros(4,fdv.ntp);

    %Manage trivial case of vb=0: no signal, or vb=1: all vessel
%     if (vb<eps) 
%         %All return variables have been initialized to zero, above.
%         return
%     elseif (vb>(1-eps))
%         % If segmented acquisition with NFlips>NTP, where VIF has been
%         % measured or interpolated to each excitation
%         % Combine segments: ** NOTE HERE ASSUMES EQUAL WEIGHTING **
%         for ii=1:fdv.ntp 
%             Mxyiv(:,ii)=mean(MxyivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
%             Mziv(:,ii)=mean(MivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
%         end
%         return
%     end
    
    
    %Calculate nontrivial result
    %
    % Diff Eq: y'(t) = Ay(t) + ff(t)
    % Sol'n:   y(t) = exp(A*t)*integral(0,t: exp(-A*tau)*ff(tau) dtau)
    %
    %Here:
    % [Pee;Lee;Pic;Lic] = [Pee0;Lee0;Pic0;Lic0]*exp(At) + ...
    %                                    kpl*integral(0,t:exp(A*(t-tau))*[Piv(tau);Liv(tau);0's;0's]*dtau)
    % A is 2x2
    a11 = -(kvedve+kecpdve+(R1Pyr));
    a13 = kecpdve;
    a22 = -(kvedve+kecldve+(R1Lac));
    a24 = kecldve;
    a31 = kecpdvc;
    a33 = -(kecpdvc+x.kpl+(R1Pyr));
    a34 = x.klp;
    a42 = kecldvc;
    a43 = x.kpl;
    a44 = -(kecldvc+x.klp+(R1Lac));

    A = [a11 0 a13 0; 0 a22 0 a24; a31 0 a33 a34; 0 a42 a43 a44];
    % Diagonlize to permit matrix integral.  
    [P,D]=eig(A);
    % D is diagonlized matrix (2x2)
    % P is matrix of eigenvectors 
    %   --> A*P=P*D
    %   --> A=P*D\P
    dD=diag(D); % get diagonal values (1x2)
    
    %Calculate signal evolution each TR:
    
    MzevSegIC=[x.Pe0; x.Le0; x.Pi0; x.Li0]; %Set Initial Conditions
    for ii=1:fdv.ntp
        MxyevSeg=zeros(4,fdv.NSeg);
        MzevSeg =zeros(4,fdv.NSeg);
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
                MzevSeg2a=exp(dD*TR).*(-b./dD).*(exp(-dD*TR)-1);
                MzevSeg2b=exp(dD*TR).*m.*(((-TR./dD)-(1./dD./dD)).*exp(-dD*TR)+(1./dD./dD));
                
                %Total signal at end of TR equals IC for next TR:
                MzevSegIC = P*(MzevSeg1 + kvedve*(MzevSeg2a+MzevSeg2b));
            end
        end
        if (fdv.NSeg>1) %manage signal from segmented acquisitions
            %sizemxyev=size(mean(MxyevSeg))
            Mxyev(:,ii)=mean(MxyevSeg,2);
            Mzev(:,ii) =mean(MzevSeg,2);
            Mxyiv(:,ii)=mean(MxyivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
            Mziv(:,ii) =mean(MzivSeg(:,(ii-1)*fdv.NSeg+1:ii*fdv.NSeg),2);
        else
            Mxyev(:,ii)=MxyevSeg;
            Mzev(:,ii) =MzevSeg;
            Mxyiv(:,ii)=MxyivSeg(:,ii);
            Mziv(:,ii) =MzivSeg(:,ii);
        end            
    end

if fdv.verbose
    figure(99)
    plot(fdv.taxis(1:fdv.NSeg:end),Mxyev(1,:),'g-',...
        fdv.taxis(1:fdv.NSeg:end),Mxyev(3,:),'g:',...
        fdv.taxis(1:fdv.NSeg:end),Mxyev(2,:),'b-',...
        fdv.taxis(1:fdv.NSeg:end),Mxyev(4,:),'b:')
    legend('PyrEV','LacEV','PyrIC','LacIC')
end    