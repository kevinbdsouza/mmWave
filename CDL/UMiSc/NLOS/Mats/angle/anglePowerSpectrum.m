clc
clear 

%parameters
numUsers  = 10;
nF = 1;

for nU = 2:2
    
AoDSpec = load('AoDUMiNlosF');
AoDSpec = AoDSpec.AoDUMiNlosF;
AoASpec = load('AoAUMiNlosF');
AoASpec = AoASpec.AoAUMiNlosF;

%average over all the rays 
AoDSpec = squeeze(mean(AoDSpec,3));
AoASpec = squeeze(mean(AoASpec,3));

AoDperUser = AoDSpec(:,:,nU,nF);
totalPathsAod = size(AoDperUser,2);
plotAod(AoDperUser(1,:),AoDperUser(2,:),AoDperUser(3,:),totalPathsAod,nU);

AoAperUser = AoASpec(:,:,nU,nF);
totalPathsAoa = size(AoAperUser,2);
plotAoa(AoAperUser(1,:),AoAperUser(2,:),AoAperUser(3,:),totalPathsAoa,nU);

end 

function [] = plotAoa(Uaoa,Vaoa,Waoa,totalPathsAoa,nU)
Xaoa = zeros(1,totalPathsAoa);
Yaoa = zeros(1,totalPathsAoa);
Zaoa = zeros(1,totalPathsAoa);

h1 = figure;
quiver3(Xaoa,Yaoa,Zaoa,Uaoa,Vaoa,Waoa,'LineWidth',1,'ShowArrowHead','off')
set(gca,'fontsize',12)
%title(['AoA spectrum for User ',num2str(nU)])
xlabel('x')
ylabel('y')
zlabel('z')
view(-80,30)
%saveas(h1,sprintf('AoAUser%d.png',nU));
end 

function [] = plotAod(Uaod,Vaod,Waod,totalPathsAod,nU)

Xaod = zeros(1,totalPathsAod);
Yaod = zeros(1,totalPathsAod);
Zaod = zeros(1,totalPathsAod);

h2 = figure;
quiver3(Xaod,Yaod,Zaod,Uaod,Vaod,Waod,'LineWidth',1,'ShowArrowHead','off')
title(['AoD spectrum for User ',num2str(nU)])
xlabel('x')
ylabel('y')
zlabel('z')
view(20,30)
%saveas(h2,sprintf('AoDUser%d.png',nU));
end 
