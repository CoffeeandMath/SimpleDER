cd ..
cd ObjectOrientedConstrained
RunCodeEndLoad
cd ..
cd AnalyticalEndforceDef
Runner

plot([pi^2/4 pi^2/4], [-1,1])
xlim([0 20])
ylim([-1 1])

xlabel('Applied Force')
ylabel('End Tip Location')

legend('x computational','y computational', 'x analytical' , 'y analytical' ,'Euler Buckling load')
