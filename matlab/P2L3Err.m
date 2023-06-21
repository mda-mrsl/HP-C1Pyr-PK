function [resid] = P2L3Err(parms,fdv)

[EV,IV,vols] = P2L3(parms,fdv);
tot=(IV(1:2,:)*vols(1))+(EV(1:2,:)*vols(2))+(EV(3:4,:)*vols(3));
resid=(fdv.data-tot);
%Normalize to max of measured data. May need to re-weight...
resid(1,:)=resid(1,:)/max(fdv.data(1,:));
resid(2,:)=resid(2,:)/max(fdv.data(2,:));

resid=reshape(resid,1,prod(size(fdv.data))); 