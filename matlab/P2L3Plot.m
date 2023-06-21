function fh=P2L3Plot(parms,fdv)

%Unpack variables:
for ii=1:length(parms)
    x.(fdv.fitvars{ii})=parms(ii);
end
for ii=1:length(fdv.knowns)
    x.(fdv.knowns{ii})=fdv.knownvals(ii);
end

NSeg=fdv.NSeg;

[EV,IV,vols,Mzev,Mziv] = P2L3(parms,fdv);
tot=(IV(1:2,:)*vols(1))+(EV(1:2,:)*vols(2))+(EV(3:4,:)*vols(3));
%Mz =(Mziv*vols(1))+(Mzev(1:2,:)*vols(2))+(Mzev(3:4,:)*vols(3));
resid=P2L3Err(parms,fdv);
resn=sqrt(sum(resid.^2));

fh=figure;
% First plot Mz
%yyaxis left 
subplot(2,1,1)
plot(fdv.taxis(NSeg:NSeg:end),Mziv(1,:)*vols(1),'g-',...
    fdv.taxis(NSeg:NSeg:end),Mzev(1,:)*vols(2),'g-.',...
    fdv.taxis(NSeg:NSeg:end),Mzev(3,:)*vols(2),'g:',...
    fdv.taxis(NSeg:NSeg:end),Mzev(2,:)*vols(3),'b-.',...
    fdv.taxis(NSeg:NSeg:end),Mzev(4,:)*vols(3),'b:',...
    'linewidth',2)
legend('Pyr_{iv}','Pyr_{ee}','Pyr_{ic}','Lac_{ee}','Lac_{ic}','Location','Best')
ylabel('Mz')
title(sprintf('%s: L2Norm = %10.4e',fdv.Name,resn))
%yyaxis right
subplot(2,1,2)
plot(fdv.taxis(NSeg:NSeg:end),tot(1,:),'g-',...
    fdv.taxis(NSeg:NSeg:end),tot(2,:),'b-',...
    fdv.taxis(NSeg:NSeg:end),fdv.data(1,:),'gx',...
    fdv.taxis(NSeg:NSeg:end),fdv.data(2,:),'bx',...
    'linewidth',2,'MarkerSize',8)
title(sprintf('k_p_l=%5.3f, k_v_e=%5.3f, v_b=%6.3e',x.kpl,x.kve,x.vb))
ylabel('Mxy')
xlabel('Time (s)')
legend('Pyr','Lac','Location','Best')
