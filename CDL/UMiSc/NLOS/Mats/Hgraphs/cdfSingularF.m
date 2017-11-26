clc
clear 

numUsers = 10;

%------------------------6GHz--------------------------
Huser = load('HPerF6.mat');  
Huser = Huser.HPerF6;
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
axis([0 5 0 1])
title('Empirical CDF of Singular Values')
xlabel('Singular Value')
ylabel('Empirical CDF')
for timeStamp = 1:1
singularAvgMat1 = singularTime1(:,:,timeStamp);
singularAvg1 = mean(singularAvgMat1,1);
[f1,x1] = ecdf(singularAvg1);
end 


%--------------30GHz-------------------
Huser = load('HPerF30.mat');  
Huser = Huser.HPerF30;
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


for timeStamp = 1:1
singularAvgMat1 = singularTime1(:,:,timeStamp);
singularAvg1 = mean(singularAvgMat1,1);
[f2,x2] = ecdf(singularAvg1);
end 

%---------------------------60GHz-----------------------------
Huser = load('HPerF60.mat');  
Huser = Huser.HPerF60;
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


for timeStamp = 1:1
singularAvgMat1 = singularTime1(:,:,timeStamp);
singularAvg1 = mean(singularAvgMat1,1);
[f3,x3] = ecdf(singularAvg1);
end 

%----------------------70GHz-------------------------------
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

for timeStamp = 1:1
singularAvgMat1 = singularTime1(:,:,timeStamp);
singularAvg1 = mean(singularAvgMat1,1);
[f4,x4] = ecdf(singularAvg1);
end 

plot(x1,f1)
plot(x2,f2)
plot(x3,f3)
plot(x4,f4)


hold off
