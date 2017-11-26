clc
clear 

%parameters
numUsers = 10; 

for nU = 1:numUsers

    %AOA matrix 
    aoaFile = sprintf('AOALobePowerSpectrum%d.mat', nU);
    aoaSpecStruct = importdata(aoaFile);
    fieldsAoa = fieldnames(aoaSpecStruct);
    numLobesAoa = length(fieldsAoa);
    
    for temp = 1:numLobesAoa 
        if isempty(aoaSpecStruct.(fieldsAoa{temp})) 
            aoaSpecStruct = rmfield(aoaSpecStruct,fieldsAoa{temp});
        end 
    end 
    
    fieldsAoa = fieldnames(aoaSpecStruct);
    numLobesAoa = length(fieldsAoa);
    aoaSpec = [];
    for nL = 1:numLobesAoa
    aoaTemp = aoaSpecStruct.(fieldsAoa{nL});
    aoaSpec = cat(1,aoaSpec,aoaTemp);
    end 
    mpcAoa(aoaSpec,nU);
    
     %AOD matrix
    aodFile = sprintf('AODLobePowerSpectrum%d.mat', nU);
    aodSpecStruct = importdata(aodFile);
    fieldsAod = fieldnames(aodSpecStruct);
    numLobesAod = length(fieldsAod);
    
    for temp = 1:numLobesAod 
        if isempty(aodSpecStruct.(fieldsAod{temp})) 
            aodSpecStruct = rmfield(aodSpecStruct,fieldsAod{temp});
        end 
    end 
    
    fieldsAod = fieldnames(aodSpecStruct);
    numLobesAod = length(fieldsAod);
    aodSpec = [];
    for nL = 1:numLobesAod
    aodTemp = aodSpecStruct.(fieldsAod{nL});
    aodSpec = cat(1,aodSpec,aodTemp);
    end 
    mpcAod(aodSpec,nU);

end 

function [] = mpcAoa(aoaSpec,nU)
Uaoa = [];
Vaoa = [];
Waoa = [];
totalPathsAoa = 0;
    for nPath = 1:size(aoaSpec,1)
    aoaPathTemp = aoaSpec(nPath,:);
    power = aoaPathTemp(2);
    aoaAz = aoaPathTemp(4) * (pi/180);
    aoaEle = aoaPathTemp(5) * (pi/180);
    power = abs(10*log10(power) + 30);
    uAoa = power*sin(aoaEle)*cos(aoaAz);
    vAoa = power*sin(aoaEle)*sin(aoaAz);
    wAoa = power*cos(aoaEle);
    Uaoa = [Uaoa uAoa];
    Vaoa = [Vaoa vAoa]; %#ok<*AGROW>
    Waoa = [Waoa wAoa];
    totalPathsAoa = totalPathsAoa + 1;    
    end 
    plotAoa(Uaoa,Vaoa,Waoa,totalPathsAoa,nU);
end 

function [] = mpcAod(aodSpec,nU)
Uaod = [];
Vaod = [];
Waod = [];
totalPathsAod = 0;
    for nPath = 1:size(aodSpec,1)
    aodPathTemp = aodSpec(nPath,:);
    power = aodPathTemp(2);
    aodAz = aodPathTemp(4) * (pi/180);
    aodEle = aodPathTemp(5) * (pi/180);
    power = abs(10*log10(power) + 30);
    uAod = power*sin(aodEle)*cos(aodAz);
    vAod = power*sin(aodEle)*sin(aodAz);
    wAod = power*cos(aodEle);
    Uaod = [Uaod uAod];
    Vaod = [Vaod vAod];
    Waod = [Waod wAod];
    totalPathsAod = totalPathsAod + 1;    
    end
    plotAod(Uaod,Vaod,Waod,totalPathsAod,nU);
end 


function [] = plotAoa(Uaoa,Vaoa,Waoa,totalPathsAoa,nU)
Xaoa = zeros(1,totalPathsAoa);
Yaoa = zeros(1,totalPathsAoa);
Zaoa = zeros(1,totalPathsAoa);

h1 = figure;
quiver3(Xaoa,Yaoa,Zaoa,Uaoa,Vaoa,Waoa,'LineWidth',1,'ShowArrowHead','off')
title(['AoA spectrum for User ',num2str(nU)])
xlabel('x')
ylabel('y')
zlabel('z')
hold on
view(-80,30)
saveas(h1,sprintf('AoAUser%d.png',nU));
hold off
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
hold on
view(80,30)
saveas(h2,sprintf('AoDUser%d.png',nU));
hold off

end 
