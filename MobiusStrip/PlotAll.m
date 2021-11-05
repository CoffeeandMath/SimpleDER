
formatSpec = '%f';

etaname = 'eta.txt';
k1name = 'k1.txt';
k3name = 'k3.txt';

wvals = cell(5,1);
wvals{1} = 0.1;
wvals{2} = 0.2;
wvals{3} = 0.5;
wvals{4} = 0.8;
wvals{5} = 1;

colorlist = cell(size(wvals));
colorlist{1} = 'r';
colorlist{2} = 'g';
colorlist{3} = 'b';
colorlist{4} = 'k';
colorlist{5} = 'c';

etavals = cell(size(wvals));
k1vals = cell(size(wvals));
k3vals = cell(size(wvals));

svals = cell(size(wvals));



for i = 1:size(wvals,1)
    cd(['w_' num2str(wvals{i})])
    
    etaID = fopen(etaname,'r');
    k1ID = fopen(k1name,'r');
    k3ID = fopen(k3name,'r');
    
    eta = fscanf(etaID,formatSpec);
    k1 = fscanf(k1ID,formatSpec);
    k3 = fscanf(k3ID,formatSpec);
    
    Nn = size(eta,1);
    
    fclose(etaID);
    fclose(k1ID);
    fclose(k3ID);
    
    
    eta = circshift(eta(1:(end-1)),floor(Nn/2));
    k1 = circshift([k1],floor((Nn-1)/2));
    k3 = circshift([k3],floor((Nn-1)/2));
    
    s = linspace(0,1,size(k1,1));
    
    etavals{i} = eta;
    % Note that the 2pi normalization is due to Starostin using a beam of
    % length 2 pi
    k1vals{i} = k1*Nn/(2*pi);
    k3vals{i} = k3*Nn/(2*pi);  
    svals{i} = s;
    
    
    cd ..
end

wvalsname = cell(size(wvals));
for i = 1:size(wvals,1)
    wvalsname{i} = num2str(wvals{i});
end

figure()
hold all
for i = 1:size(wvals,1)
    
    plot(svals{i},-1.0*k1vals{i},colorlist{i})
    
end
legend(wvalsname)
title('\kappa_1')

figure()
hold all
for i = 1:size(wvals,1)
    
    plot(svals{i},k3vals{i},colorlist{i})
    
end
legend(wvalsname)
title('\kappa_3')
ylim([-3 4])

figure()
hold all
for i = 1:size(wvals,1)
    
    plot(svals{i},etavals{i},colorlist{i})
    
end
legend(wvalsname)
title('\eta')
