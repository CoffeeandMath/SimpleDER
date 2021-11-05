cd ..
cd ObjectOrientedConstrained
RunCodeGravLoad
cd ..
cd AnalyticalGravDef
Runner

plot([7.83735 7.83735], [-1,1])
xlim([0 150])
ylim([-1 1])

xlabel('Applied Force')
ylabel('End Tip Location')

legend('x computational','y computational', 'x analytical' , 'y analytical' ,'Euler Buckling load')
