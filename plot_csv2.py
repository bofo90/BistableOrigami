import matplotlib as matl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#import pylab as P
import csv
#import glob
#import ntpath
#import sys

plt.close('all')

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
    
def cm2inch(value):
    return value/2.54

def NiceGraph3D(axes, nameX, nameY, nameZ, mincoord = [np.NaN, np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN, np.NaN],
                divisions = [np.NaN, np.NaN, np.NaN], buffer = [0.0, 0.0, 0.0]):
    gray = '0.2'
    matl.rcParams.update({'font.size': 15})

    if ~np.isnan(mincoord[0]) and ~np.isnan(maxcoord[0]):
        axes.set_xlim3d([mincoord[0]-buffer[0], maxcoord[0]+buffer[0]])
        if ~np.isnan(divisions[0]):
            axes.set_xticks(np.linspace(mincoord[0],maxcoord[0],divisions[0]))
    axes.set_xlabel(nameX)
    axes.xaxis.labelpad = 20
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim3d([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if ~np.isnan(divisions[1]):
            axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
    axes.set_ylabel(nameY)
    axes.yaxis.labelpad = 20
    if ~np.isnan(mincoord[2]) and ~np.isnan(maxcoord[2]):
        axes.set_zlim3d([mincoord[2]-buffer[2], maxcoord[2]+buffer[2]])
        if ~np.isnan(divisions[2]):
            axes.set_zticks(np.linspace(mincoord[2],maxcoord[2],divisions[2]))
    axes.set_zlabel(nameZ)
    axes.zaxis.labelpad = 15
    
    axes.w_xaxis.set_pane_color((0,0,0,0))
    axes.w_yaxis.set_pane_color((0,0,0,0))
    axes.w_zaxis.set_pane_color((0,0,0,0))
    
    axes.w_xaxis.line.set_color(gray)
    axes.w_yaxis.line.set_color(gray)
    axes.w_zaxis.line.set_color(gray)
    
    axes.w_xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray)
    axes.w_yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray)
    axes.w_zaxis.label.set_color(gray)
    axes.tick_params(axis='z', colors=gray)    
    return

def NiceGraph2D(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
                buffer = [0.0, 0.0, 0.0]):
    gray = '0.2'
    matl.rcParams.update({'font.size': 15})

    if ~np.isnan(mincoord[0]) and ~np.isnan(maxcoord[0]):
        axes.set_xlim([mincoord[0]-buffer[0], maxcoord[0]+buffer[0]])
        if isinstance(divisions[0], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[0]).any():
                axes.set_xticks(divisions[0])
        else:
            if ~np.isnan(divisions[0]):
                axes.set_xticks(np.linspace(mincoord[0],maxcoord[0],divisions[0]))
                            
    axes.set_xlabel(nameX)
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if ~isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
    axes.set_ylabel(nameY)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    return

#######################################################################################################################
##################### Initializing data from files 
#######################################################################################################################
hingeNum = np.empty(0, dtype = int)
hingeName = np.empty(0)
eEdgeFol = np.empty(0)
eEdgeRel = np.empty(0)
eFaceFol = np.empty(0)
eFaceRel = np.empty(0)
eHingeFol = np.empty(0)
eHingeRel = np.empty(0)
eTAngleFol = np.empty(0)
eTAngleRel = np.empty(0)
eHinIntFol = np.empty(0)
eHinIntRel = np.empty(0)

exflFol = np.empty(0, dtype = int)
exflRel = np.empty(0, dtype = int)

CMxFol = np.empty(0)
CMxRel = np.empty(0)
CMyFol = np.empty(0)
CMyRel = np.empty(0)
CMzFol = np.empty(0)
CMzRel = np.empty(0)
RadFol = np.empty(0)
RadRel = np.empty(0)
StdFol = np.empty(0)
StdRel = np.empty(0)
MaxStrFol = np.empty(0)
MaxStrRel = np.empty(0)
MinStrFol = np.empty(0)
MinStrRel = np.empty(0)

folder_name = "Results/cube/sqp/energy/"

file_name1 = "EnergyData.csv"
file_name2 = "Hinges.csv"
file_name3 = "PosStad.csv"

#######################################################################################################################
##################### Reading Files
#######################################################################################################################
with open(folder_name + file_name1,'r',newline='') as ffit: 
    readers = csv.reader(ffit,delimiter=',');
    for row in readers:
        if row:
            if isfloat(row[0].encode('ascii','ignore')):
                hingeNum = np.append(hingeNum, np.int(row[0]));
                eEdgeFol = np.append(eEdgeFol, np.float(row[1]));
                eEdgeRel = np.append(eEdgeRel, np.float(row[2]));
                eFaceFol = np.append(eFaceFol, np.float(row[3]));
                eFaceRel = np.append(eFaceRel, np.float(row[4]));
                eHingeFol = np.append(eHingeFol, np.float(row[5]));
                eHingeRel = np.append(eHingeRel, np.float(row[6]));
                eTAngleFol = np.append(eTAngleFol, np.float(row[7]));
                eTAngleRel = np.append(eTAngleRel, np.float(row[8]));
                eHinIntFol = np.append(eHinIntFol, np.float(row[9]))
                eHinIntRel = np.append(eHinIntRel, np.float(row[10]))
                exflFol = np.append(exflFol, np.int(row[11]))
                exflRel = np.append(exflRel, np.int(row[12]))
                
with open(folder_name + file_name2,'r',newline='') as ffit: 
    readers = csv.reader(ffit,delimiter=',');
    for row in readers:
        if row:
            if isfloat(row[0].encode('ascii','ignore')):
                hingeName = np.append(hingeName, row[1]);

with open(folder_name + file_name3,'r',newline='') as ffit: 
    readers = csv.reader(ffit,delimiter=',');
    for row in readers:
        if row:
            if isfloat(row[0].encode('ascii','ignore')):
                CMxFol = np.append(CMxFol, np.float(row[2]));
                CMxRel = np.append(CMxRel, np.float(row[3]));
                CMyFol = np.append(CMyFol, np.float(row[4]));
                CMyRel = np.append(CMyRel, np.float(row[5]));
                CMzFol = np.append(CMzFol, np.float(row[6]));
                CMzRel = np.append(CMzRel, np.float(row[7]));
                RadFol = np.append(RadFol, np.float(row[8]));
                RadRel = np.append(RadRel, np.float(row[9]));
                StdFol = np.append(StdFol, np.float(row[10]));
                StdRel = np.append(StdRel, np.float(row[11]));
                MaxStrFol = np.append(MaxStrFol, np.float(row[12]));
                MaxStrRel = np.append(MaxStrRel, np.float(row[13]));
                MinStrFol = np.append(MinStrFol, np.float(row[14]));
                MinStrRel = np.append(MinStrRel, np.float(row[15]));
                
ffit.close()

#######################################################################################################################
##################### Analize the data
#######################################################################################################################

stepsHinge = int(len(hingeNum)/hingeNum[-1])
totalflags = 6
internalHinges = 12 ###### Number of internal hinges

tolHinge = 0.0001
tolEdge = 0.01

############################################ get the number of actuated hinges for each hinge-set
actuatedHinges = np.zeros(hingeNum[-1], dtype = int)
hingeCount = np.zeros(internalHinges-1)   
for hinge in np.arange(hingeNum[-1]):
    actuatedHinges[hinge] = len(hingeName[hinge].split())
    hingeCount[actuatedHinges[hinge]-1] += 1
############################################ order the hinge-set according to the number of actuated hinges
orderedHinges = np.argsort(actuatedHinges)

############################################ count the different flags and Mask the hinge-sets that didnt converged
flagCountFol = np.zeros((len(hingeCount), totalflags))
flagCountRel = np.zeros((len(hingeCount), totalflags))
hingesMask = np.arange(hingeNum[-1], dtype = float)
notConvHinges = np.empty((0,3), dtype = int)
for i in np.arange(len(hingeNum)):
    if exflFol[i] != 1 or exflRel[i] != 1:
        flagCountFol[actuatedHinges[hingeNum[i]-1]-1,exflFol[i]+3]  += 1
        flagCountRel[actuatedHinges[hingeNum[i]-1]-1,exflRel[i]+3]  += 1
        notConvHinges = np.append(notConvHinges,np.array([[hingeNum[i]-1, exflFol[i], exflRel[i]]]), axis = 0)
        hingesMask[hingeNum[i]-1] = np.NaN
notConverged = sum(np.isnan(hingesMask))
converged = hingeNum[-1] - notConverged
############################################ normalize the flag counts
for i in np.arange(totalflags):
    flagCountFol[:,i] = flagCountFol[:,i]/hingeCount
    flagCountRel[:,i] = flagCountRel[:,i]/hingeCount
    
############################################ ordering the last energies on the release according to the number of actuated hinges
lasteHingeRel = eHingeRel[stepsHinge-1::stepsHinge]
lasteEdgeRel = eEdgeRel[stepsHinge-1::stepsHinge]
lasteHingeRel = lasteHingeRel[orderedHinges]
lasteEdgeRel = lasteEdgeRel[orderedHinges]
energies = np.append(lasteHingeRel, lasteEdgeRel)
############################################ going through all hinges in order and see all the neighbours in energy
#nonConvStates = 0
differentEnergies = np.empty((0,2), dtype = int)
for hinge in np.arange(hingeNum[-1]):
    if ~np.isnan(hingesMask[orderedHinges[hinge]]):
        hingeeHinge = lasteHingeRel[hinge]
        hingeeEdge = lasteEdgeRel[hinge]
        sameEnergy = np.where(np.logical_and(np.logical_and(energies[:hingeNum[-1]]>=hingeeHinge-tolHinge,
                                                            energies[:hingeNum[-1]]<=hingeeHinge+tolHinge),
                                             np.logical_and(energies[hingeNum[-1]:]>=hingeeEdge-tolEdge, 
                                                            energies[hingeNum[-1]:]<=hingeeEdge+tolEdge)))[0]
        ############################################ see if at least one of the neighbours converge
        i = 0
        if len(notConvHinges) > 0:
            while len(np.where(notConvHinges[:,0] == orderedHinges[sameEnergy[i]])[0]) != 0:
                i = i + 1
                if i >= len(sameEnergy):
                    i = -1
                    break
        ############################################ save the new state or not according to the convergance
        if i != -1:
            findit = np.where(differentEnergies[:,0] == orderedHinges[sameEnergy[i]])[0]
            if len(findit) == 0:
                differentEnergies = np.append(differentEnergies,np.array([[orderedHinges[sameEnergy[i]], len(sameEnergy)]]), axis = 0)
        else:
#            nonConvStates += len(sameEnergy)
            print('State from non convergent hinges')
#differentEnergies = np.array(differentEnergies)


#######################################################################################################################
##################### Ploting the result
#######################################################################################################################
fig1 = plt.figure(0,figsize=(cm2inch(35), cm2inch(20)))
ax1 = plt.subplot(111)#
NiceGraph2D(ax1, 'Edge Energy',  'Hinge Energy')#, [min(eEdgeRel[stepsHinge-1::stepsHinge]), np.NaN],[0.0014,  np.NaN] )
#ax3 = plt.subplot(122)
#NiceGraph2D(ax3, 'Edge Energy',  'Hinge Energy', [0.00025, 0.012],[0.01375,  0.029] )

fig2 = plt.figure(1,figsize=(cm2inch(35), cm2inch(20)))
ax2 = plt.subplot(111)
NiceGraph2D(ax2, 'Average Radius', 'StDev of Radius')#, [min(eEdgeRel[stepsHinge-1::stepsHinge]), np.NaN],[max(eEdgeRel[stepsHinge-1::stepsHinge]), np.NaN], buffer = [0.00001, 0.0])

fig3 = plt.figure(2,figsize=(cm2inch(35), cm2inch(20)))
ax4 = plt.subplot(111, projection='3d')
NiceGraph3D(ax4, 'CenterMass X', 'CenterMass Y', 'CenterMass Z')

fig4 = plt.figure(3,figsize=(cm2inch(35), cm2inch(20)))
ax5 = plt.subplot(111)
NiceGraph2D(ax5, 'Hinge-Set Number', 'Internal Hinge Energy', [np.nan, min(eHinIntRel)], [np.nan, max(eHinIntRel)], buffer = [0, 0.0004])

fig5 = plt.figure(4,figsize=(cm2inch(35), cm2inch(20)))
ax6 = plt.subplot(121)  
ax7 = plt.subplot(122)   
maxPerc = np.amax([np.nanmax(flagCountFol),np.nanmax(flagCountRel)])
NiceGraph2D(ax6, '# of actuated Hinges', 'percentage # of flags', [0.5, 0], [internalHinges-0.5, 1+0.01], [np.arange(len(hingeCount))+1, np.nan])       
NiceGraph2D(ax7, '# of actuated Hinges', 'percentage # of flags', [0.5, 0], [internalHinges-0.5, 1+0.01], [np.arange(len(hingeCount))+1, np.nan])  

fig6 = plt.figure(5,figsize=(cm2inch(35), cm2inch(20)))
ax8 = plt.subplot(121)  
ax9 = plt.subplot(122) 
NiceGraph2D(ax8, 'Hinge-Set Number', 'Max final streching')
NiceGraph2D(ax9, 'Hinge-Set Number', 'Min final streching')            
            
width = 0.3
separation = 0
colors = cm.Set2(np.linspace(0, 1, totalflags))
for i, c in zip(reversed(np.arange(totalflags)), reversed(colors)):
    if i == 4 or i == 0 or i == 2:  ###### block to appear flags 1, -3 and -1 respectively
        continue
    ax6.bar(np.arange(len(hingeCount))+separation, flagCountFol[:,i], width, color=c, label = i-3)
    ax7.bar(np.arange(len(hingeCount))+separation, flagCountRel[:,i], width, color=c)
    separation += width

ax6.legend(loc = 2, fontsize =15, framealpha = 0.5, edgecolor = 'inherit', fancybox = False)
          
colors = cm.rainbow(np.linspace(0, 1, hingeNum[-1]))
for hinge, c in zip(np.arange(hingeNum[-1]),colors):
    if ~np.isnan(hingesMask[hinge]):
        c = '#36648B'
    else:
        c = '#FE9128'
    if len(np.where(differentEnergies[:,0] == hinge)[0]) != 0:
        ax1.scatter(eEdgeRel[stepsHinge*hinge+stepsHinge-1], eHingeRel[stepsHinge*hinge+stepsHinge-1], c = c, label = hingeName[hinge])
#        ax3.scatter(eEdgeRel[stepsHinge*hinge+stepsHinge-1], eHingeRel[stepsHinge*hinge+stepsHinge-1], c = c)
        ax2.scatter(RadRel[stepsHinge*hinge+stepsHinge-1], StdRel[stepsHinge*hinge+stepsHinge-1], c = c, label = hingeName[hinge])
        ax4.scatter(CMxRel[stepsHinge*hinge+stepsHinge-1], CMyRel[stepsHinge*hinge+stepsHinge-1], CMzRel[stepsHinge*hinge+stepsHinge-1], c = c)
        ax5.scatter(hingeNum[stepsHinge*hinge+stepsHinge-1], eHinIntRel[stepsHinge*hinge+stepsHinge-1], c = c)
        ax8.scatter(hingeNum[stepsHinge*hinge+stepsHinge-1], MaxStrRel[stepsHinge*hinge+stepsHinge-1], c = c)
        ax9.scatter(hingeNum[stepsHinge*hinge+stepsHinge-1], abs(MinStrRel[stepsHinge*hinge+stepsHinge-1]), c = c)
#        ax1.plot(eEdgeRel[stepsHinge*hinge:stepsHinge*(hinge+1)],  eHingeRel[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = c)
#        ax3.plot(eEdgeRel[stepsHinge*hinge:stepsHinge*(hinge+1)],  eHingeRel[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = c)
#        ax2.plot(RadRel[stepsHinge*hinge:stepsHinge*(hinge+1)],  StdRel[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = c)
    if len(np.where(differentEnergies[:,0] == hinge)[0]) != 0:
        ax1.annotate(hingeName[hinge], xy=(eEdgeRel[stepsHinge*hinge+stepsHinge-1], eHingeRel[stepsHinge*hinge+stepsHinge-1]), 
                      xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
                      arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
        
        



fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()
fig6.tight_layout()
fig1.show()
fig2.show()
fig3.show()
fig4.show()
fig5.show()
fig6.show()
fig1.savefig(folder_name + file_name1[:-4]+'.png', transparent = True)
#fig2.savefig(folder_name + file_name3[:-4]+'.png', transparent = True)
#fig3.savefig(folder_name + 'CenterOfMass.png', transparent = True)
fig5.savefig(folder_name + 'Flags.png', transparent = True)
#fig6.savefig(folder_name + 'MaxMinStretch.png', transparent = True)