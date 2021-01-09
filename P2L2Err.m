    function [resid] = P2L2Err(parms,fdv)

[EV,IV,vb] = P2L2(parms,fdv);

tot=(IV*vb)+EV*(1-vb);

resid=(fdv.data-tot);

%Normalize to max of measured data. May need to adjust weights...
resid(1,:)=resid(1,:)/max(fdv.data(1,:));
resid(2,:)=resid(2,:)/max(fdv.data(2,:));

resid=reshape(resid,1,prod(size(fdv.data))); 