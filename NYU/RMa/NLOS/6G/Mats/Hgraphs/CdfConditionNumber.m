clc
clear 

%numUsers = 10;
nUser = 4;
Huser = load('HuserUMiNlos6.mat');
Huser = Huser.Huser;
nSub = size(Huser,3);
%user 4 and 5 condition numbers above 80 (pretty high)
HperUser = Huser(:,:,:,nUser);
condNum = zeros(1,nSub);

for nK = 1:nSub
    HperSubk = HperUser(:,:,nK);
    condTemp = cond(HperSubk); 
    condNum(nK) = condTemp;
end

condNumdB = mag2db(condNum);
[f,x] = ecdf(condNumdB);
plot(x,f)
axis([0 120 0 1])
title('Empirical CDF of Condition number')
xlabel('Condition number (dB)')
ylabel('Empirical CDF')
