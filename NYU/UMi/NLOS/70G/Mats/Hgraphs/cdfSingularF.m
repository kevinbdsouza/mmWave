clc
clear 

numUsers = 10;

%-------------6GHz-------------
Huser = load('HUMiNlos6NYU.mat');  
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
axis([0 8 0 1])
title('Empirical CDF of Singular Values')
xlabel('Singular Value')
ylabel('Empirical CDF')
singularAvg1 = mean(singularAvgMat1,1);
[f1,x1] = ecdf(singularAvg1);

%----------------------30GHz--------------------------------
Huser = load('HUMiNlos30NYU.mat');  
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

singularAvg1 = mean(singularAvgMat1,1);
[f2,x2] = ecdf(singularAvg1);


%-------------------------60GHz------------------
Huser = load('HUMiNlos60NYU.mat');  
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


singularAvg1 = mean(singularAvgMat1,1);
[f3,x3] = ecdf(singularAvg1);

%---------------------70GHz-----------------------------
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


singularAvg1 = mean(singularAvgMat1,1);
[f4,x4] = ecdf(singularAvg1);untitled


plot(x1,f1)
plot(x2,f2)
plot(x3,f3)
plot(x4,f4)
set(gca,'fontsize',12)
legend('6GHz','30GHz','60GHz','70GHz')

hold off



