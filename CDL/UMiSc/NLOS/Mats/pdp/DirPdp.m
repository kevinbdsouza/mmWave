
pdp = load('pdpUMiNlosF');
pdp = pdp.pdpUMiNlosF;

numUsers = 10;
nF = 4;

for nU = 1:numUsers
pdpPerUser = pdp(:,:,nU,nF);

h1 = figure;

h=stem(pdpPerUser(:,1)*10^9, pdpPerUser(:,2), '.');
set(h,'BaseValue',-40);
title(['Omnidirectional Pdp for User ',num2str(nU)])
xlabel('Delay (ns)')
ylabel('Received power (dB)')
grid on
saveas(h1,sprintf('OmniDirPdp%d.png',nU));

end 