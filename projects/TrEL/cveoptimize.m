% Calculate CVE and optimize for Whittaker Smooth
function [zout, cveval, lambdaval, pks, locs, optimumlambda] = cveoptimize(datainput,minlambdaexp, maxlambdaexp,weightvector)
lambdaval = logspace(minlambdaexp, maxlambdaexp);
counter = 1;
tempcveval = zeros(1,50);
tempz = zeros(1,length(datainput));
if nargin < 4
    for i=lambdaval
        [~, cve,  ~] = whitsm(datainput, i,2);
        tempcveval(counter) = cve;
        counter = counter + 1;
    end
elseif nargin == 4
    for i = lambdaval
        [~, cve, ~] = whitsmw(datainput, weightvector, i,2);
        tempcveval(counter) = cve;
        counter = counter + 1;
    end
end
deriv = gradient(tempcveval);
deriv2 = gradient(deriv);
%subplot(3,1,1)
%semilogx(lambdaval, tempcveval)
%subplot(3,1,2)
%semilogx(lambdaval,deriv)
%hold on
%semilogx(lambdaval,deriv2)
%hold off
[temppks, templocs] = findpeaks(deriv2);
%for i = templocs
%    xline(lambdaval(i))
%end
firstpeak = templocs(1);
if tempcveval(firstpeak)>min(tempcveval)
    optimumlambda = lambdaval(firstpeak);
else
    optimumlambda = min(tempcveval);
end
if nargin < 4
    [zout, ~, ~] = whitsm(datainput, optimumlambda, 2);
elseif nargin == 4
    [zout, ~, ~] = whitsmw(datainput, weightvector, optimumlambda, 2);
end
%subplot(3,1,3)
%plot(datainput,'r*')
%hold on
%plot(zout,'b')
cveval = tempcveval;
pks = temppks;
locs = templocs;
end