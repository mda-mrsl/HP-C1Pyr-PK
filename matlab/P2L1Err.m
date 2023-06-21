function [resid] = P2L1Err(parms,fdv)

[Mxy] = P2L1(parms,fdv);
resid=(fdv.data(2,:)-Mxy(2,:));

%Normalize to max of measured data.
resid=resid/max(fdv.data(2,:));
