clc
clear 

numUsers = 10;
Huser = load('Huser.mat');  
Huser = Huser.Huser;
nSub = size(Huser,3);

singularAvgMat1 = [];
singularAvgMat2 = [];
for nUser = 1:numUsers
HperUser = Huser(:,:,:,nUser);
singular1 = zeros(1,nSub);
singular2 = zeros(1,nSub);

for nK = 1:nSub
    HperSubk = HperUser(:,:,nK);
    %HperSubk = HperSubk.';
    S = svd(HperSubk); 
    singular1(nK) = S(1);
    singular2(nK) = S(2);
end 
singularAvgMat1 = [singularAvgMat1;singular1]; %#ok<*AGROW>
singularAvgMat2 = [singularAvgMat2;singular2];
end



figure(); hold on
axis([0 10 0 1])
title('Empirical CDF of Singular Values')
xlabel('Singular Value')
ylabel('Empirical CDF')
%condNumAvgMat = condNumAvgMat([1:3,6:10],:);  %removes very poor users 
singularAvg1 = mean(singularAvgMat1,1);
singularAvg2 = mean(singularAvgMat2,1);
%singularAvgdB1 = mag2db(singularAvg1);
%singularAvgdB2 = mag2db(singularAvg2);
[f1,x1] = ecdf(singularAvg1);
[f2,x2] = ecdf(singularAvg2);
plot(x1,f1)
plot(x2,f2)

