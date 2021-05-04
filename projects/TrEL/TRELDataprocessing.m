% Batch Processing
% Assumes data is stored in a Cell named 'workingcell' with the offset
% voltages, and applied voltages stored in row 2 col 7,8

offsetvoltages = workingcell{2,7};
appliedvoltages = workingcell{2,8};
tempmatrix = zeros(6,5);
tempmatrix(1,:) = appliedvoltages;
expnum = 5; %Sets the upper limit for traversing the cell
delayrow = 6; %Delay Time data should be in this row
outputcell = workingcell;
for i = 1:expnum %Main Loop (Traverses cell)
    for j = 1:length(appliedvoltages)
        tempdata = workingcell{6,i+1}(j,:);
        [curve, h, k] = fit(offsetvoltages',tempdata','power2');
        tempcoeff = coeffvalues(curve);
        tempmatrix(2,j) = tempcoeff(1);% a value
        tempmatrix(3,j) = tempcoeff(2);% b value
        tempmatrix(4,j) = tempcoeff(3);% c value
        tempmatrix(5,j) = h.rsquare;
        %a*x^b + c = y
        % logb(-c/a) = xint
        tempmatrix(6,j) = log(-tempcoeff(3)/tempcoeff(1))/log(tempcoeff(2));
    end
    outputcell{8,i+1} = tempmatrix;
    tempmatrix = zeros(6,5);
    tempmatrix(1,:) = appliedvoltages;
end