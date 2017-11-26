clc
clear 

nRuns = 2;
frequency = [6e9,30e9,60e9,70e9];
pdpUMiNlosF = [];
AoDUMiNlosF = [];
AoAUMiNlosF = [];
HUMiNlosF = [];

for fc = 1:length(frequency)

pdpRun = [];
AoDRun = [];
AoARun = [];
HRun = []; 
    
for run = 1:nRuns 
    
pd = makedist('Uniform','lower',10,'upper',200/3);
d2D = random(pd,1,2);
d2D = ceil(d2D);
pdpUsers = [];
AoDUsers = [];
AoAUsers = [];
HUsers = [];
    
for i = 1:length(d2D)
[UMiScNlosChan,UMiScNlosChanH,txWaveform,Param] = UMiScNlos(d2D(i),frequency(fc));
[PDP,AoDSpec,AoASpec] = UMiScNlosChan(txWaveform);
[HperUser,~,~] = UMiScNlosChanH(txWaveform);
pdpUsers = cat(3,pdpUsers,PDP(:,1:2));
AoDUsers = cat(4,AoDUsers,AoDSpec);
AoAUsers = cat(4,AoAUsers,AoASpec);
chanMat = squeeze(sum(HperUser,2));
HUsers = cat(5,HUsers,chanMat);
end 

pdpRun = cat(4,pdpRun,pdpUsers);
AoDRun = cat(5,AoDRun,AoDUsers);
AoARun = cat(5,AoARun,AoAUsers);
HRun = cat(6,HRun,HUsers); 

end 

pdpPerF = mean(pdpRun,4);
AoDPerF = mean(AoDRun,5);
AoAPerF = mean(AoARun,5);

if fc == 1
HPerF6  = mean(HRun,6);
elseif fc == 2
HPerF30  = mean(HRun,6);   
elseif fc == 3
HPerF60  = mean(HRun,6);  
elseif fc == 4
HPerF70  = mean(HRun,6); 
end 

pdpUMiNlosF = cat(4,pdpUMiNlosF,pdpPerF);
AoDUMiNlosF = cat(5,AoDUMiNlosF,AoDPerF);
AoAUMiNlosF = cat(5,AoAUMiNlosF,AoAPerF);

end 

%coherenceTime = 1/(2 * Param.maxDopShift);
