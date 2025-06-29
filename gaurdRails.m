function [newPset] = gaurdRails(Pset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here: index of specific variable/param +
%   lower bound + upper bound
global NAtot Atot AMPKtot
Rails=[2 0 AMPKtot; 11 0 NAtot-Pset(13)-Pset(12); 12 0  NAtot-Pset(11)-Pset(13); 13 0  NAtot-Pset(11)-Pset(12); 6 0 Atot];%6:AMP; 2:pAMPK; 11:NAM
newPset=Pset;
for i=1:length(Rails(:,1))
    param=Rails(i,1);
    if Pset(param)<Rails(i,2)
        newPset(param)=Rails(i,2);
    elseif Pset(param)>Rails(i,3)
        newPset(param)=Rails(i,3);
    end
end
end

