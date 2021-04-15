% Clear all variables
clear
% Prompt User to pick TREL files
files = uipickfiles;
% Set the timescale. Should be consistent for these files
timescale = 83.335;
datacell=cell(15,length(files)+1);
lambdaval = logspace(-2,8);
cveval = zeros(1,length(lambdaval));
t0ind = 120;
counter = 1;
difforder = 2; %Change this value to change the smoothing
%headers
datacell{1,1} = 'filename - 1';
datacell{2,1} = 'rawdata - 2';
datacell{3,1} = 'index - 3';
datacell{4,1} = 't95 - 4';
datacell{5,1} = 't50 - 5';
datacell{6,1} = 'Smoothed ROI data - 6';
datacell{7,1} = 'Derivative data - 7';
datacell{8,1} = 'Slope - 8';
datacell{9,1} = 'Delay Time - 9';
datacell{10,1} = 'Linear Time - 10';
datacell{11,1} = 'Lambda Value - 11';
datacell{12,1} = 'Smoothing Order - 12';
datacell{13,1} = 'Voltages - 13';
datacell{14,1} = 'OffsetVoltages - 14';
datacell{15,1} = 'Fall Time - 15';
for i = 2:length(files)+1
    [a,name,b]=fileparts(files{i-1});
    datacell{1,i}=name;
    datacell{2,i}=readmatrix(files{i-1});
    datacell{2,i}(:,1)=(datacell{2,i}(:,1)*timescale)-10000.2;
    datacell{2,i}(:,2)=datacell{2,i}(:,2)*(1/10000); %normalized, and zereod data
    datacell{3,i}=i; %index of each file
    t99ind = find(datacell{2,i}(:,2)>0.99,1,'first');
    t95ind = find(datacell{2,i}(:,2)>0.95,1,'first');
    t50ind = find(datacell{2,i}(:,2)>0.5,1,'first');
    fallind = find(datacell{2,i}(:,1)>=100000,1,'first');
    fallt05ind = fallind + find(datacell{2,i}(fallind:end,2)<=0.05,1,'first');
    fallwindow = abs(fallind - fallt05ind);
    falltime = datacell{2,i}(fallt05ind,1)-datacell{2,i}(fallind,1);
    t99 = datacell{2,i}(t99ind,1);
    t95 = datacell{2,i}(t95ind,1);
    t50 = datacell{2,i}(t50ind,1);
    datacell{4,i}=t95;
    datacell{5,i}=t50;
    windowsize = abs(t95ind-t0ind);
    Xorig = datacell{2,i}(t0ind:t95ind,1);
    Yorig = datacell{2,i}(t0ind:t95ind,2);
    XYorig = datacell{2,i}(t0ind:t95ind,1:2);
    if windowsize<101 %Check to see if the data window is large enough 
        for o = 1:6       
        modwindow(o) = (2^o)*windowsize - (2^o) - 1;
        end
        %Find how many iterations of resample need to be performed
        iterations = find(modwindow>100,1,'first');
        if windowsize<51
            windowexpansion = ceil((51-windowsize)/2);
            XYorig = datacell{2,i}((t0ind-windowexpansion):(t95ind+windowexpansion),1:2);
            iterations = 1;
        end
        [resampled, weightvector] = resampling(XYorig,iterations);
        Xresample = resampled(:,1); 
        [Ysmooth, cveval, lambdaval,optpk,optloc,optimumlambda] = cveoptimize(resampled(:,2),-2,8,weightvector);
        temp = zeros(length(Xresample),2);
        temp(:,1) = Xresample;
        temp(:,2) = Ysmooth;
    elseif windowsize>101
        [Ysmooth, cveval, lambdaval,optpk,optloc,optimumlambda] = cveoptimize(Yorig,-2,8);
        temp = zeros(length(Xorig),2);
        temp(:,1)=Xorig;
        temp(:,2)=Ysmooth;
    end
    datacell{6,i} = temp;
    rawderiv = gradient(Yorig)./gradient(Xorig);
    if windowsize<101
        smoothderiv = gradient(Ysmooth)./gradient(Xresample);
    elseif windowsize>101
        smoothderiv = gradient(Ysmooth)./gradient(Xorig);
    end
    datacell{7,i} = smoothderiv;
    [pks, locs] = findpeaks(smoothderiv);
    maxpkind = find(pks==max(pks));
    slope = pks(maxpkind);
    datacell{8,i}=slope;
    if windowsize<101
        lineartime = Xresample(locs(maxpkind));
        linearvalue = Ysmooth(locs(maxpkind));
    elseif windowsize>101
        lineartime = Xorig(locs(maxpkind));
        linearvalue = Ysmooth(locs(maxpkind));
    end
    yint = linearvalue - slope * lineartime;
    tdelay = -yint/slope;
    datacell{9,i} = tdelay;
    datacell{10,i} = lineartime;
    datacell{11,i} = optimumlambda;
    datacell{12,i} = difforder;
    datacell{15,i} = falltime;
    if windowsize<101
        yline = slope*Xresample + yint;
    elseif windowsize>101
        yline = slope*Xorig + yint;
    end
    subplot(3,2,1)
    plot(datacell{2,i}(:,1),datacell{2,i}(:,2))
    xline(t95,'-',num2str(t95))
    xline(t50,'-',num2str(t50))
    axis([-2000 110000 0 1])
    title('Raw Data')
    xlabel('time (ns)')
    ylabel('Normalized Intensity (a.u.)')
    subplot(3,2,2)
    semilogx(lambdaval,cveval)
    for p = optloc
        xline(lambdaval(p),'-',num2str(cveval(p)))
    end
    title('Cross Validation Error')
    xlabel('lambda (a.u.)')
    ylabel('CVE (a.u.)')
    subplot(3,2,3)
    plot(Xorig,Yorig,'*b')
    hold on
    if windowsize<101
        plot(Xresample,Ysmooth)
    elseif windowsize>101
        plot(Xorig, Ysmooth)
    end
    hold off
    title('Smoothed & Original Data')
    xlabel('time (ns)')
    ylabel('Normalized Intensity (a.u.)')
    axis([0-timescale,t95+timescale,0,1])
    subplot(3,2,4)
    plot(Xorig,rawderiv);
    hold on
    if windowsize<101
        plot(Xresample,smoothderiv);
    elseif windowsize>101
        plot(Xorig,smoothderiv);
    end
    hold off
    xline(lineartime,'-',{num2str(lineartime),num2str(slope)})
    title('Derivative of Raw and Smoothed Data')
    xlabel('time (ns)')
    ylabel('Slope (I/ns)')
    axis tight
    subplot(3,2,5)
    plot(Xorig,Yorig)
    hold on
    if windowsize<101
        plot(Xresample,yline)
    elseif windowsize>101
        plot(Xorig,yline)
    end
    xline(tdelay,'-',num2str(tdelay))
    xline(lineartime,'-',num2str(lineartime))
    xline(t50,'-',num2str(t50))
    xline(t95,'-',num2str(t95))
    axis([0-83.335, t95+83.335, 0, 1,])
    hold off
    title('Raw Data + Linear Region')
    xlabel('time (ns)')
    ylabel('Normalized Intensity (a.u.)')
    subplot(3,2,6)
    if windowsize<101
        plot(Xresample,Ysmooth,Xresample,yline)
    elseif windowsize>101
        plot(Xorig,Ysmooth,Xorig,yline)
    end
    axis([0-83.335, t95+83.335, 0, 1])
    title('Smoothed Data + Linear Region')
    xlabel('time (ns)')
    ylabel('Normalized Intensity (a.u.)')
    set(gcf,'position',[100, 40, 1080, 1300])
    savename = strcat(datacell{1,i},'.png');
    saveas(gcf,savename)
    %pause
    clf
    counter = 1;
    clear resampled
end

voltages = input("Voltages [V1,V2,...]: ");
offsetvoltages = input("Offset Voltages [V1,V2,...]: ");
datacell{13,2} = voltages;
datacell{14,2} = offsetvoltages;

writecell(datacell);

% Additional Plotting
rows = length(voltages);
col = length(offsetvoltages);
datacelldim = size(datacell);
numcol = datacelldim(2);
fomrow = [4,5,8,9,10];
for i = 1:rows
    for j = 1:col
        t95mat(i,j) = datacell{4,((j+1)+(col*(i-1)))};
    end
end
for i = 1:rows
    for j = 1:col
        t50mat(i,j) = datacell{5,((j+1)+(col*(i-1)))};
    end
end
for i = 1:rows
    for j = 1:col
        slopemat(i,j) = datacell{8,((j+1)+(col*(i-1)))};
    end
end
for i = 1:rows
    for j = 1:col
        delaymat(i,j) = datacell{9,((j+1)+(col*(i-1)))};
    end
end
for i = 1:rows
    for j = 1:col
        tlinearmat(i,j) = datacell{10,((j+1)+(col*(i-1)))};
    end
end
for i = 1:rows
    for j = 1:col
        falltimemat(i,j) = datacell{15,((j+1)+(col*(i-1)))};
    end
end
plotcell = {t95mat,t50mat,slopemat,delaymat,tlinearmat,falltimemat};
titlecell = {'t95', 't50', 'Slope', 'Delay Time', 'Linear Time','Fall Time'};
legendlist = cell(1,col);
for n = 1:col
    legendlist{n} = num2str(offsetvoltages(n));
end
for k = 1:length(plotcell)
    subplot(3,2,k)
    hold on
    for l = 1:col
        plot(voltages,plotcell{k}(:,l))
    end
    title(titlecell{k})
    legend(legendlist)
    xlabel('Applied Voltage')
    ylabel('Time (ns)')
    hold off
end
set(gcf,'position',[100, 40, 1080, 1300])
pause
saveas(gcf,'extracted-data.png')
writecell(plotcell,'extracted-data')
        