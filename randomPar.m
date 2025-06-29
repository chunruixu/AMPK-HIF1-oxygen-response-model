function [paramSet] = randomPar(paramSet,paramTarg,stdv)
%[paramSet] = samplePar(paramSet,paramTarg,stdv)
%Returns a new parameter set with variation around a normal distribution
%for parameters specified by paramTarg
for i=1:length(paramTarg)
    paramSet(paramTarg(i))=max(paramSet(paramTarg(i))+paramSet(paramTarg(i))*normrnd(0,stdv/100),0);
    paramSet=gaurdRails(paramSet);  %restrictions for boundary of paraSet
end