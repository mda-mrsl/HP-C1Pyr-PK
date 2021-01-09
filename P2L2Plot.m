function fh=P2L2Plot(parms,fdv)

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

[EV,IV,vb,Mzev,Mziv] = P2L2(parms,fdv);
tot=(IV*vb)+EV*(1-vb);
Mz =(Mziv*vb)+Mzev*(1-vb);
resid=P2L2Err(parms,fdv);
resn=sqrt(sum(resid.^2));

fh=gcf;
% First plot Mz
%yyaxis left 
subplot(2,1,1)
plot(fdv.taxis(NSeg:NSeg:end),Mz(Pyr,:),'g-',...
    fdv.taxis(NSeg:NSeg:end),Mz(Lac,:),'b-',...
    'linewidth',2)
legend('Pyr','Lac','Location','Best')
ylabel('Mz')
title(sprintf('%s: L2Norm = %10.4e',fdv.Name,resn))
%yyaxis right
subplot(2,1,2)
plot(fdv.taxis(NSeg:NSeg:end),tot(Pyr,:),'g-',...
    fdv.taxis(NSeg:NSeg:end),tot(Lac,:),'b-',...
    fdv.taxis(NSeg:NSeg:end),fdv.data(Pyr,:),'gx',...
    fdv.taxis(NSeg:NSeg:end),fdv.data(Lac,:),'bx',...
    'linewidth',2,'MarkerSize',8)
title(sprintf('k_p_l=%5.3f, k_v_e=%5.3f, v_b=%6.3e',x.kpl,x.kve,x.vb))
ylabel('Mxy')
xlabel('Time (s)')
legend('Pyr','Lac','Location','Best')
