clc
clear 

numUsers = 10;
Huser = load('Huser.mat');
Huser = Huser.Huser;
nSub = size(Huser,3);
condNumAvgMat = [];
for nUser = 1:numUsers
 %best = 3
HperUser = Huser(:,:,:,nUser);
condNum = zeros(1,nSub);

for nK = 1:nSub
    HperSubk = HperUser(:,:,nK);
    condTemp = cond(HperSubk); 
    condNum(nK) = condTemp;
end 
condNumAvgMat = [condNumAvgMat;condNum]; %#ok<*AGROW>
end

%condNumAvgMat = condNumAvgMat([1:3,6:10],:);  %removes very poor users   
condNumAvg = mean(condNumAvgMat,1);

condNumAvgdB = mag2db(condNumAvg);
[f,x] = ecdf(condNumAvgdB);
plot(x,f)
axis([0 120 0 1])
title('Empirical CDF of Condition number')
xlabel('Condition number (dB)')
ylabel('Empirical CDF')
