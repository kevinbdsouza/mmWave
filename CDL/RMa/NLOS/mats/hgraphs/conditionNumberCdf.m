clc
clear 

numUsers = 10;
Huser = load('HPerF70.mat');  
Huser = Huser.HPerF70;
nSub = size(Huser,4);
   
condNumAvgMat = [];
for nUser = 1:numUsers
HperUser = squeeze(Huser(1,:,:,:,nUser));
condNum = zeros(1,nSub);

for nK = 1:nSub
    HperSubk = HperUser(:,:,nK);
    HperSubk = HperSubk.';
    condTemp = cond(HperSubk); 
    condNum(nK) = condTemp;
end 
condNumAvgMat = [condNumAvgMat;condNum]; %#ok<*AGROW>
end

%condNumAvgMat = condNumAvgMat([1:3,6:10],:);  %removes very poor users 
condNumAvg = mean(condNumAvgMat,1);
condNumAvgdB = mag2db(condNumAvg);
[f,x] = ecdf(condNumAvgdB);
figure(); 
plot(x,f)
axis([0 40 0 1])
title('Empirical CDF of Condition number')
xlabel('Condition number (dB)')
ylabel('Empirical CDF')

