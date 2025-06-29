function score = findSquares2(t,y,data)
yvalues=zeros(1,length(data(:,1)));

for i=1:length(data(:,1))
    ID = findClosest(t, data(i,1));
    yvalues(i) = y(ID);
end

score = sum((yvalues'-data(:,2)).^2)./length(data(:,1));
