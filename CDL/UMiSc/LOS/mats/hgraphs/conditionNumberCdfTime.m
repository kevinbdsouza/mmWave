clc
clear 

numUsers = 10;
Huser = load('HPerF30.mat');  
Huser = Huser.HPerF30;
nSub = size(Huser,4);
condNumTime = []; 

for timeStamp = 1:size(Huser,1)
    
condNumAvgMat = [];
for nUser = 1:numUsers
HperUser = squeeze(Huser(timeStamp,:,:,:,nUser));
condNum = zeros(1,nSub);

for nK = 1:nSub
    HperSubk = HperUser(:,:,nK);
    HperSubk = HperSubk.';
    condTemp = cond(HperSubk); 
    condNum(nK) = condTemp;
end 
condNumAvgMat = [condNumAvgMat;condNum]; %#ok<*AGROW>
end
condNumTime = cat(3,condNumTime,condNumAvgMat); 
end 



figure(); hold on
axis([0 40 0 1])
title('Empirical CDF of Condition number')
xlabel('Condition number (dB)')
ylabel('Empirical CDF')
for timeStamp = 1:size(Huser,1)
%condNumAvgMat = condNumAvgMat([1:3,6:10],:);  %removes very poor users 
condNumAvgMat = condNumTime(:,:,timeStamp);
condNumAvg = mean(condNumAvgMat,1);
condNumAvgdB = mag2db(condNumAvg);
[f,x] = ecdf(condNumAvgdB);
plot(x,f)
end 
hold off
