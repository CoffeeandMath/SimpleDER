%% Plotting phi vs y0 norm
close all
a = 0.2;

bvals = [0.2,0.4,0.6,0.8,1,1.5,2,8];

N = max(size(bvals));

figure()
for i = 1:N
    load(['a_' num2str(a) '_b_' num2str(bvals(i)) '_phivsy0.mat']);
    plot(phivals*180/pi,-y0*pi)
    hold all
end
%plot(phivals*180/pi,-(1-2/pi)*phivals)

legnames = cell(N,1);
for i = 1:N
    legnames{i} = ['Aspect Ratio b/a = ' num2str(bvals(i)/a)];
end
xticks([0 60 120 180])
xticklabels({'0','60','120','180'})
yticks([-0.4 0 .4])
legend(legnames)