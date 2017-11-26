clc
clear 

nUser = 2;    
Huser = load('HuserUMiLos70.mat');
Huser = Huser.Huser;
HperUser = Huser(:,:,:,nUser);
nSub = size(HperUser,3);
condNum = zeros(1,nSub);

for nK = 1:nSub
    HperSubk = HperUser(:,:,nK);
    condTemp = cond(HperSubk); 
    condNum(nK) = condTemp;
end 

condNumdB = mag2db(condNum);
[f,x] = ecdf(condNumdB);
plot(x,f)
axis([0 50 0 1])
title('Empirical CDF of Condition number')
xlabel('Condition number (dB)')
ylabel('Empirical CDF')
