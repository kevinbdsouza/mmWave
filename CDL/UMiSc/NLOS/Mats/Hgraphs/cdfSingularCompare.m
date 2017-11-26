clc
clear 

numUsers = 10;

%CDL
Huser = load('HPerF70.mat');  
Huser = Huser.HPerF70;
nSub = size(Huser,4);
singularTime1 = []; 
singularTime2 = [];

for timeStamp = 1:size(Huser,1)
    
singularAvgMat1 = [];
singularAvgMat2 = [];
for nUser = 1:numUsers
HperUser = squeeze(Huser(timeStamp,:,:,:,nUser));
singular1 = zeros(1,nSub);
singular2 = zeros(1,nSub);

for nK = 1:nSub
    HperSubk = HperUser(:,:,nK);
    HperSubk = HperSubk.';
    S = svd(HperSubk); 
    singular1(nK) = S(1);
    singular2(nK) = S(2);
end 
singularAvgMat1 = [singularAvgMat1;singular1]; %#ok<*AGROW>
singularAvgMat2 = [singularAvgMat2;singular2];
end
singularTime1 = cat(3,singularTime1,singularAvgMat1); 
singularTime2 = cat(3,singularTime2,singularAvgMat2);
end


figure(); hold on
axis([0 2.5 0 1])
title('Empirical CDF of Singular Values')
xlabel('Singular Value')
ylabel('Empirical CDF')
for timeStamp = 1:1
%condNumAvgMat = condNumAvgMat([1:3,6:10],:);  %removes very poor users 
singularAvgMat1 = singularTime1(:,:,timeStamp);
singularAvg1 = mean(singularAvgMat1,1);
singularAvgMat2 = singularTime2(:,:,timeStamp);
singularAvg2 = mean(singularAvgMat2,1);
%singularAvgdB1 = mag2db(singularAvg1);
%singularAvgdB2 = mag2db(singularAvg2);
[f1,x1] = ecdf(singularAvg1);
[f2,x2] = ecdf(singularAvg2);
plot(x1,f1)
plot(x2,f2)
end


%----------------NYU---------------------
Huser = load('HUMiNlos70NYU.mat');  
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


%condNumAvgMat = condNumAvgMat([1:3,6:10],:);  %removes very poor users 
singularAvg1 = mean(singularAvgMat1,1);
singularAvg2 = mean(singularAvgMat2,1);
%singularAvgdB1 = mag2db(singularAvg1);
%singularAvgdB2 = mag2db(singularAvg2);
[f1,x1] = ecdf(singularAvg1);
[f2,x2] = ecdf(singularAvg2);
plot(x1,f1)
plot(x2,f2)
hold off
set(gca,'fontsize',12)
legend('3GPP Larger','3GPP Smaller','NYU Larger','NYU Smaller')
