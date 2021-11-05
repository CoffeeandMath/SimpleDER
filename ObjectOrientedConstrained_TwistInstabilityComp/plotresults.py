import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

#plt.style.use('seaborn-whitegrid')
plt.style.use('default')
fig, axs = plt.subplots(nrows=2, ncols=2,sharex=True)
DataThin = pd.read_csv('DataThin.csv')
DataThick = pd.read_csv('DataThick.csv')
SanoDataThinExp = pd.read_csv('SanoThinDatay0Exp.csv')
SanoDataThinSim = pd.read_csv('SanoThinDatay0Sim.csv')
SanoDataThickExp = pd.read_csv('SanoThickDatay0Exp.csv')
SanoDataThickSim = pd.read_csv('SanoThickDatay0Sim.csv')

#Defining plot colors and markers
y0DERcolor = '#ad2e2f'
y0SIMcolor = '#b7b072'
y0EXPcolor = '#b7b072'

SIMmarker = 's'
EXPmarker = '^'

eigDERcolor = '#ad2e2f'
DERLabel = 'Discrete elastic rod simulations'
SANOSIMLabel = 'Sano et al. simulations'
SANOEXPLabel = 'Sano et al. experiments'
SnapthroughLabel = 'Snap-through'
SnapshotLabel = 'Snapshots of configuration'

#ax1.plot(EndLoadAnalytical['Force'], EndLoadAnalytical['X1'], label='Analytical Solution',color=AnalyticalLinecolor,linestyle=AnalyticalLinestyle,linewidth = LineWidth)
offset1 = -1
snapsplit = 636
axs[0,0].plot(180.*DataThin['phi'][:snapsplit]/np.pi,DataThin['y0'][:snapsplit],color=y0DERcolor,label=DERLabel,clip_on=False)
axs[0,0].plot(180.*DataThin['phi'][(snapsplit+1):offset1]/np.pi,DataThin['y0'][(snapsplit+1):offset1],color=y0DERcolor,clip_on=False)
axs[1,0].plot(180.*DataThin['phi'][:snapsplit]/np.pi,350.*DataThin['minEig'][:snapsplit], color = eigDERcolor,clip_on=False)
axs[1,0].plot(180.*DataThin['phi'][(snapsplit+1):offset1]/np.pi,350.*DataThin['minEig'][(snapsplit+1):offset1], color = eigDERcolor,clip_on=False)

offset2 = 1;
axs[0,1].plot(180.*DataThick['phi']/np.pi,DataThick['y0'],color=y0DERcolor,clip_on=False)
axs[1,1].plot(180.*DataThick['phi']/np.pi,350.*DataThick['minEig'], color = eigDERcolor,clip_on=False)

#Plotting the image points
indices1 = [0,249,749,998]
indices2 = [0,249,749,992]

axs[0,0].scatter(180.*DataThin['phi'][indices1]/(np.pi), DataThin['y0'][indices1],color = y0DERcolor,marker = 'o',clip_on=False, label = SnapshotLabel)
axs[1,0].scatter(180.*DataThin['phi'][indices1]/(np.pi), 350.*DataThin['minEig'][indices1],color = y0DERcolor,marker = 'o',clip_on=False)

axs[0,1].scatter(180.*DataThick['phi'][indices2]/(np.pi), DataThick['y0'][indices2],color = y0DERcolor,marker = 'o',clip_on=False)
axs[1,1].scatter(180.*DataThick['phi'][indices2]/(np.pi), 350.*DataThick['minEig'][indices2],color = y0DERcolor,marker = 'o',clip_on=False)


SanoDataXthinResc1 = SanoDataThinExp['xUnsc'] - SanoDataThinExp['xUnsc'][0]
SanoDataXthinResc2 = SanoDataXthinResc1*180./SanoDataXthinResc1[len(SanoDataXthinResc1)-1]
SanoDataYthinResc1 = SanoDataThinExp['yUnsc'] - SanoDataThinExp['yUnsc'][0]
SanoDataXthinRescSim1 = SanoDataThinSim['xUnsc'] - SanoDataThinSim['xUnsc'][0]
SanoDataXthinRescSim2 = SanoDataXthinRescSim1*180./SanoDataXthinRescSim1[len(SanoDataXthinRescSim1)-1]
SanoDataYthinRescSim1 = SanoDataThinSim['yUnsc'] - SanoDataThinSim['yUnsc'][0]
y0Data = 138.04
y04Data = 115.66
SanoDataYthinResc2 = 0.4*SanoDataYthinResc1/(y0Data - y04Data)
SanoDataYthinRescSim2 = 0.4*SanoDataYthinRescSim1/(y0Data - y04Data)
axs[0,0].scatter(SanoDataXthinResc2,SanoDataYthinResc2,marker = EXPmarker, facecolors='none', edgecolors=y0EXPcolor,label = SANOEXPLabel)
axs[0,0].scatter(SanoDataXthinRescSim2,SanoDataYthinRescSim2,marker = SIMmarker, facecolors='none', edgecolors=y0SIMcolor, label = SANOSIMLabel)



SanoDataXThickRescSim1 = SanoDataThickSim['xUnsc'] - SanoDataThickSim['xUnsc'][0]
SanoDataXThickRescSim2 = SanoDataXThickRescSim1*180./SanoDataXThickRescSim1[len(SanoDataXthinRescSim1)-1]
SanoDataYThickRescSim1 = SanoDataThickSim['yUnsc'] - SanoDataThickSim['yUnsc'][0]
y0Data = 154.72
y08Data = 109.73
SanoDataYThickRescSim2 = 0.8*SanoDataYThickRescSim1/(y0Data - y08Data)

SanoDataYThickRescExp1 = SanoDataThickExp['yUnsc'] - SanoDataThickSim['yUnsc'][0]
SanoDataYThickRescExp2 = 0.8*SanoDataYThickRescExp1/(y0Data - y08Data)
#axs[0,1].plot(SanoDataXthinResc2,SanoDataYthinResc2,linestyle=None,marker = 'o')

#print(len(SanoDataYThickRescExp2))
axs[0,1].scatter(SanoDataXThickRescSim2[:13],SanoDataYThickRescExp2,marker = EXPmarker, facecolors='none', edgecolors=y0EXPcolor,clip_on=False)
axs[0,1].scatter(SanoDataXThickRescSim2,SanoDataYThickRescSim2,marker = SIMmarker, facecolors='none', edgecolors=y0SIMcolor,clip_on=False)

instInd = 635
axs[0,0].axvline(180.*(DataThin['phi'][instInd]+DataThin['phi'][instInd+1])/(2.0*np.pi), color = 'green', linestyle=':', linewidth = 3, label = SnapthroughLabel,clip_on=False)
axs[1,0].axvline(180.*(DataThin['phi'][instInd]+DataThin['phi'][instInd+1])/(2.0*np.pi), color = 'green', linestyle=':', linewidth = 3,clip_on=False)






axs[0,0].set_xlim([0, 180])
#axs[1,1].set_ylim([0, 1E-7])
plt.xticks(np.arange(0, 180+1, 60.0))
axs[1,0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
axs[1,1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
fig.set_size_inches(8,6.)
#plt.subplots_adjust(wspace=0.8)
fig.legend(loc="lower center")
plt.show()
