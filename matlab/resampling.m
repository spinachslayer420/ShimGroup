% Resample time series data
function [output,weight] = resampling(datainput, iteration)
data = datainput;
dim = size(data);
len = dim(1);
rownum = (2^iteration)*len-(2^iteration - 1);
tempweight = zeros(rownum,1);
for i = 1:iteration
dim = size(data);
len = dim(1);
temp = zeros((2*len)-1,2);
    for j = 1:len
        temp((2*j-1),1)=data(j,1); %Store original times in every other row
        temp((2*j-1),2)=data(j,2); %Store original data in every other row
        if j ~= len
        temp((2*j),1) = (data(j,1)+data(j+1,1))/2;
        temp((2*j),2) = 0;
        end
    end
data = temp;
end
spacing = 2^iteration;
for i = 1:spacing:rownum
    tempweight(i) = 1;
end
output = temp;
weight = tempweight;
end