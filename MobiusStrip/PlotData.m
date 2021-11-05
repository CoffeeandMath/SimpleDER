
formatSpec = '%f';

etaID = fopen('eta.txt','r');
k1ID = fopen('k1.txt','r');
k3ID = fopen('k3.txt','r');

eta = fscanf(etaID,formatSpec);
k1 = fscanf(k1ID,formatSpec);
k3 = fscanf(k3ID,formatSpec);


fclose(etaID);
fclose(k1ID);
fclose(k3ID);

Nn = size(eta,1);
eta = circshift(eta(1:(end-1)),floor(Nn/2));
k1 = circshift([0;k1],floor(Nn/2));
k3 = circshift([0;k3],floor(Nn/2));

s = linspace(0,1,size(eta,1));

figure()
plot(s,eta)
title('\eta')
figure()
subplot(2,1,1)
plot(s,-1.0*k1*Nn/(2*pi))
title('\kappa_1')
subplot(2,1,2)
plot(s,k3*Nn/(2*pi))
title('\kappa_3')