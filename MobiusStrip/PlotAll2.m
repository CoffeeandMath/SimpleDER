cd w_eta
etafiles = dir();
wvals = [];
wvalsstr = cell(0);
cntr = 0;
for i = 3:(size(etafiles))
    ntemp = etafiles(i).name;
    numstr = ntemp(3:(end-8));
    wtemp = str2double(numstr);
    if not(isempty(wtemp))
        cntr = cntr + 1;
        
        wvals = [wvals;wtemp];
        wvalsstr{cntr} = numstr;
        
    end
end


[wvals, sval] = sort(wvals);
wvalsstrtemp = wvalsstr;
for i = 1:size(wvalsstrtemp,2)
    wvalsstr{i} = wvalsstrtemp{sval(i)};
end

fprintf('Finished cataloging w values \n');

Nwvals = max(size(wvals));
fileID = fopen(['w_' wvalsstr{1} '_eta.txt'],'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
fclose(fileID);

Neta = max(size(A));
etavals = zeros(Nwvals,Neta);

for i = 1:Nwvals
    fileID = fopen(['w_' wvalsstr{i} '_eta.txt'],'r');
    formatSpec = '%f';
    A = fscanf(fileID,formatSpec);
    fclose(fileID);
    for j = 1:Neta
        etavals(i,j) = A(j);
    end
end
cd ..

fprintf('Finished cataloging eta values \n');

cd w_k1
fileID = fopen(['w_' wvalsstr{1} '_k1.txt'],'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
fclose(fileID);

Nk1 = max(size(A));
k1vals = zeros(Nwvals,Nk1);
for i = 1:Nwvals
    fileID = fopen(['w_' wvalsstr{i} '_k1.txt'],'r');
    formatSpec = '%f';
    A = fscanf(fileID,formatSpec);
    fclose(fileID);
    for j = 1:Nk1
        k1vals(i,j) = A(j);
    end
end
cd ..
fprintf('Finished cataloging k1 values \n');
cd w_k3
fileID = fopen(['w_' wvalsstr{1} '_k3.txt'],'r');
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
fclose(fileID);

Nk3 = max(size(A));
k3vals = zeros(Nwvals,Nk3);
for i = 1:Nwvals
    fileID = fopen(['w_' wvalsstr{i} '_k3.txt'],'r');
    formatSpec = '%f';
    A = fscanf(fileID,formatSpec);
    fclose(fileID);
    for j = 1:Nk3
        k3vals(i,j) = A(j);
    end
end
cd ..
fprintf('Finished cataloging k3 values \n');
%%
w = zeros(5,1);
w(1) = 0.1/pi;
w(2) = 0.2/pi;
w(3) = 0.5/pi;
w(4) = 0.8/pi;
w(5) = 1/pi;

colorlist = cell(size(w));
colorlist{1} = 'r';
colorlist{2} = 'g';
colorlist{3} = 'b';
colorlist{4} = 'k';
colorlist{5} = 'c';

ix = zeros(size(w));



for i = 1:max(size(w))
    [ d, ix(i) ] = min( abs( wvals-w(i) ) );
    d
end


figure()
hold all
for i = 1:max(size(w))
    k1loc = circshift(k1vals(ix(i),1:(end-1)),floor((Nk1-1)/2));
    svals = linspace(0,1,max(size(k1loc)));
    plot(svals,-k1loc*(Nk1+1),colorlist{i});
end
xlim([0,1])

figure()
hold all
for i = 1:max(size(w))
    k3loc = circshift(k3vals(ix(i),1:(end-1)),floor((Nk3-1)/2));
    svals = linspace(0,1,max(size(k3loc)));
    plot(svals,k3loc*(Nk1+1),colorlist{i});
end
xlim([0,1])

figure()
hold all
for i = 1:max(size(w))
    etaloc = circshift(etavals(ix(i),:),floor(Neta/2));
    svals = linspace(0,1,max(size(etaloc)));
    plot(svals,-etaloc*(Nk1+1),colorlist{i});
end
xlim([0,1])