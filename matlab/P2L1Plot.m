function fh=P2L1Plot(parms,fdv)

Pyr=1;
Lac=2;

%Unpack variables:
for ii=1:length(parms)
    x.(fdv.fitvars{ii})=parms(ii);
end
for ii=1:length(fdv.knowns)
    x.(fdv.knowns{ii})=fdv.knownvals(ii);
end

NSeg=fdv.NSeg;

[Mxy,Mtot] = P2L1(parms,fdv);
%tot --> Mxy
%Mz --> Mtot
resid=P2L1Err(parms,fdv);
resn=sqrt(sum(resid.^2));

fh=figure;
% First plot Mz
%yyaxis left 
subplot(2,1,1)
plot(fdv.taxis(NSeg:NSeg:end),Mtot(Pyr,:),'g-',...
    fdv.taxis(NSeg:NSeg:end),Mtot(Lac,:),'b-',...
    'linewidth',2)
legend('Pyr','Lac','Location','Best')
ylabel('Mz = Mtot')
title(sprintf('%s: L2Norm = %10.4e',fdv.Name,resn))
%yyaxis right
subplot(2,1,2)
plot(fdv.taxis(NSeg:NSeg:end),Mxy(Pyr,:),'g-',...
    fdv.taxis(NSeg:NSeg:end),Mxy(Lac,:),'b-',...
    fdv.taxis(NSeg:NSeg:end),fdv.data(Pyr,:),'gx',...
    fdv.taxis(NSeg:NSeg:end),fdv.data(Lac,:),'bx',...
    'linewidth',2,'MarkerSize',8)
title(sprintf('k_p_l=%5.3f',x.kpl))
ylabel('Mxy')
xlabel('Time (s)')
legend('Pyr','Lac','Location','Best')
