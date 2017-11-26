clc
clear 


dirPdpTemp = load('DirPDPInfo.mat');  
dirPdpTemp = dirPdpTemp.DirPDPInfo;

%simulation parameters 
numUsers = 10;
numSubCar = 51;
txPower = 35; %in dbm
count = ones(1,numUsers);

%separate based on users 
for i = 1:size(dirPdpTemp,1) 
    for k = 1:numUsers
    if dirPdpTemp(i,1) == k
        row = count(1,k);
        dirPdpPerUser(row,:,k) = dirPdpTemp(i,:); %#ok<*SAGROW>
        count(1,k) = count(1,k) + 1;
    end
    end 
end 

%clean matrix 
dirPdpPerUser = dirPdpPerUser(:,3:end,:);

%calculate coeff
numTx = 64;   %URA
xNumTx = sqrt(numTx);
yNumTx = sqrt(numTx);
numRx = 2;    %ULA
carFreq = 60e9;
interCarSpacing = 2e6;
dTx = 0.5;
dRx = 0.5;

Huser = [];
for nU = 1:numUsers
    
    Hsub = [];
    for nK = 1:numSubCar
    
    subCarFreq = carFreq + (nK-1)*interCarSpacing;
    dirClean = dirPdpPerUser(:,:,nU);
    dirClean( all(~dirClean,2), : ) = [];
    
    Hpath = zeros(numRx,numTx);
    for nPath = 1:size(dirClean,1)
        dirPath = dirClean(nPath,:);
        delay = dirPath(1)*1e-9;
        rxPower = dirPath(2)+30;
        phase = dirPath(3);
        azAod = dirPath(4);
        eleAod = dirPath(5);
        azAoa = dirPath(6);
        eleAoa = dirPath(7);
        pathLoss = dirPath(8)-30;
        rmsDs = dirPath(9);
        
        rxPowerWatt = 10^(rxPower/10);
        pathGainWatt = sqrt(rxPowerWatt)*2e2;
        
        %Rectangular array response vector for TX
        indTx = 1;
        for xTx = 0:(xNumTx-1)
            for yTx = 0:(yNumTx-1)
            arrayResponseTxTemp(indTx) = exp(1i*2*pi*dTx * (xTx*sin(azAod)*sin(eleAod)+yTx*cos(eleAod)) ); 
            indTx = indTx + 1;
            end 
        end 
        arrayResponseTx = arrayResponseTxTemp.';
        
        %Linear array response vector for RX
        for xRx = 0:numRx-1
            arrayResponseRxTemp(xRx+1) = exp(1i*2*pi*dRx* xRx*sin(azAoa) ); 
        end 
        arrayResponseRx = arrayResponseRxTemp.';
        
        Hpath = Hpath + (pathGainWatt * exp(1i*phase) * exp(-1i*2*pi*subCarFreq*delay) * (arrayResponseRx * arrayResponseTx')); 
    end 
    Hsub = cat(3,Hsub,Hpath);
    
    end
 
    Huser = cat(4,Huser,Hsub);
end 



