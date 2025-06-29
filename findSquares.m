function [score,scalary] = findSquares(t,y, data, scale)
% [score,scalary] = findSquares(TOUT2,YOUT2(:,13),dpNADH);
yvalues=zeros(1,length(data(:,1)));
% cellcycletime= t(end);
% data(:,1)=data(:,1).*(cellcycletime/150);
for i=1:length(data(:,1))
    leftside=y(t-data(i,1)<=0);
%     rightside=y(t-data(i,1)>=0);
    yvalues(i)=leftside(end);
end

fun= @(scale) sum(((yvalues'-scale*data(:,2))./max(yvalues)).^2./length(data(:,1)));
if nargin < 4
    scalary = fminsearch(fun,max(yvalues));
    score = fun(scalary);
else
    score= sum(((yvalues'-scale*data(:,2))/scale).^2./length(data(:,1)));
end
end

