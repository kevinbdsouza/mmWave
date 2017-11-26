Fc = [6]; %in GHz
d2D = [35:10:500/3];
nRuns = 200;

h = figure;
hold on
for f = 1:size(Fc,2)
plAvg = [];   
for nR = 1:nRuns 
plArray = [];
for nL = 1:size(d2D,2) 
[PL] = pathloss(d2D(nL),Fc(f));
plArray = [plArray PL]; %#ok<*AGROW>
end 
plAvg = [plAvg;plArray];
end 

plAvg = mean(plAvg,1);
title('Pathloss vs Distance')
xlabel('T-R separation (m)')
ylabel('Pathloss (dB)')
plot(d2D,plAvg)
set(gca,'fontsize',12)
legend('6GHz','30GHz','60GHz','70GHz')
end 
saveas(h,sprintf('Pathloss RMa LoS.png'));

