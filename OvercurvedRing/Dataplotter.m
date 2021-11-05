%% Script to plot the D vs H results
close all
clear all
Olower = load('Olistlower.txt');
Dlower = load('Dlistlower.txt');


Oupper = load('Olistupper.txt');
Dupper = load('Dlistupper.txt');

[Osortlower,lowerInd] = sort(Olower);
Dsortlower = Dlower(lowerInd);

[Osortupper,upperInd] = sort(Oupper);
Dsortupper = Dupper(upperInd);




O = [Osortlower;Osortupper];
D = [Dsortlower;Dsortupper];
D = D/max(D);


[~,minind] = min(abs(D));
for i=minind:max(size(D))
    D(i) = -D(i);
end



plot(O,D)