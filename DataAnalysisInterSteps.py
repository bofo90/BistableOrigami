import matplotlib as matl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import from_levels_and_colors
import numpy as np
import configparser
import os.path
#import pylab as P
#import csv
#import glob
#import ntpath
#import sys

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

def ReadandAnalizeFile(folder_name, plot = True, khinge = np.nan, kedge = np.nan):
#    
#plot = True
#khinge = 0.01
#kedge = 10**0.5
#folder_name = "Results/truncated tetrahedron/sqp/energy/kh0.010_kta1.000_ke3.162/"

    plt.close('all')
    #%%
    #######################################################################################################################
    ##################### Reading Files
    #######################################################################################################################
    
    file_name1 = "/EnergyData.csv"
    file_name2 = "/Hinges.csv"
    file_name3 = "/PosStad.csv"
    file_name4 = "/Angles.csv"

    dataEnergy = np.loadtxt(folder_name+file_name1,skiprows=1, delimiter = ',', unpack = True)
    hingeNum = dataEnergy[0,:].astype(int)
    eEdgeFol = dataEnergy[1,:]
    eEdgeRel = dataEnergy[2,:]
    eFaceFol = dataEnergy[3,:]
    eFaceRel = dataEnergy[4,:]
    eHingeFol = dataEnergy[5,:]
    eHingeRel = dataEnergy[6,:]
    eTAngleFol = dataEnergy[7,:]
    eTAngleRel = dataEnergy[8,:]
    eHinIntFol = dataEnergy[9,:]
    eHinIntRel = dataEnergy[10,:]
    exflFol = dataEnergy[11,:].astype(int)
    exflRel = dataEnergy[12,:].astype(int)
    
    hingeName = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [1], dtype=bytes).astype(str)    
    
    dataPosStad = np.loadtxt(folder_name+file_name3,skiprows=1, delimiter = ',', unpack = True)
    #### If there is no metadata file the order of this data is not the same as here shown
    CMxFol = dataPosStad[1,:]
    CMxRel = dataPosStad[2,:]
    CMyFol = dataPosStad[3,:]
    CMyRel = dataPosStad[4,:]
    CMzFol = dataPosStad[5,:]
    CMzRel = dataPosStad[6,:]
    RadFol = dataPosStad[7,:]
    RadRel = dataPosStad[8,:]
    StdFol = dataPosStad[9,:]
    StdRel = dataPosStad[10,:]
    MaxStrFol = dataPosStad[11,:]
    MaxStrRel = dataPosStad[12,:]
    MinStrFol = dataPosStad[13,:]
    MinStrRel = dataPosStad[14,:]
    SumIntAngFol = dataPosStad[15,:]
    SumIntAngRel = dataPosStad[16,:]
    SumExtAngFol = dataPosStad[17,:]
    SumExtAngRel = dataPosStad[18,:]
    
    dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',')
    dataAngles = np.delete(dataAngles, 0, 1)
    MaxAngles = np.max(dataAngles, axis = 1)
    MinAngles = np.min(dataAngles, axis = 1)
    
    #%%
    #######################################################################################################################
    ##################### Read Metadata File and define variables
    #######################################################################################################################
    
    metadataFile = '/metadata.txt'
    
    if os.path.isfile(folder_name+metadataFile):
        metadata = configparser.RawConfigParser()
        metadata.read(folder_name+metadataFile)

        internalHinges = int(metadata.get('extUnitCell', 'intHinges'))
        totalnumberHinges = int(metadata.get('extUnitCell', 'Hinges'))
        totalnumberEdges = int(metadata.get('extUnitCell', 'Edges'))
        khinge = float(metadata.get('options','kHinge'))
        kedge = float(metadata.get('options','kEdge'))
        kface = float(metadata.get('options','kFace'))
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 
        #### If there is no metadata file the order of the dataPosStad is not the same as here shown           
        
        
    tolHinge = 0.003
    tolEdge = 0.01
    normalized = ''   
    tolStretch = 0.3 #max precentage of allowed stretch
    digitPi = 4 # digits of pi for rounding to see if angles go beyond it and "not converge"
    tolAngleSS = 0.087 # equivalent to 5 degrees
        
    stepsHinge = int(len(hingeNum)/len(hingeName))
    totHingeNum = len(hingeName)
    totalflags = 8
    #%%
    #######################################################################################################################
    ##################### Analize the data
    #######################################################################################################################


    ############################################ Modify data to have a normalized energy (or just difference in angle/length)
    if ~np.isnan(khinge):
        eHingeFol = np.sqrt(eHingeFol*2/khinge/totalnumberHinges)
        eHingeRel = np.sqrt(eHingeRel*2/khinge/totalnumberHinges)
        eHinIntFol = np.sqrt(eHinIntFol*2/khinge/internalHinges)
        eHinIntRel = np.sqrt(eHinIntRel*2/khinge/internalHinges)
        tolHinge = 0.008
        normalized = normalized + 'hn'
    if ~np.isnan(kedge):
        eEdgeFol = np.sqrt(eEdgeFol*2/kedge/totalnumberEdges)
        eEdgeRel = np.sqrt(eEdgeRel*2/kedge/totalnumberEdges)
        tolEdge = 0.0001
        normalized = normalized + 'en'
    
    ############################################ get the number of actuated hinges for each hinge-set
    actuatedHinges = np.zeros(totHingeNum, dtype = int)
    hingeCount = np.zeros(internalHinges)   
    for hinge in np.arange(totHingeNum):
        actuatedHinges[hinge] = len(hingeName[hinge].split())
        hingeCount[actuatedHinges[hinge]-1] += 1
    ############################################ order the hinge-set according to the number of actuated hinges
    orderedHinges = np.argsort(actuatedHinges)
    
    ############################################ count the different flags and Mask the hinge-sets that didnt converge
    flagCountFol = np.zeros((internalHinges, totalflags))
    flagCountRel = np.zeros((internalHinges, totalflags))
    hingesMask = np.arange(totHingeNum, dtype = float)
    notConvHinges = np.empty((0,3), dtype = int)
    for i in np.arange(totHingeNum):
        error = False
#        flagCountFol[actuatedHinges[i]-1,exflFol[i*stepsHinge+1]+3]  += 1
#        flagCountRel[actuatedHinges[i]-1,exflRel[i*stepsHinge+1]+3]  += 1
        if exflFol[(i+1)*stepsHinge-1] != 1:# and exflFol[(i+1)*stepsHinge-1] != 2:
            flagCountFol[actuatedHinges[i]-1,exflFol[(i+1)*stepsHinge-1]+3]  += 1
            error = True
        if exflRel[(i+1)*stepsHinge-1] != 1 and not error:# and exflRel[(i+1)*stepsHinge-1] != 2:
            flagCountRel[actuatedHinges[i]-1,exflRel[(i+1)*stepsHinge-1]+3]  += 1
            error = True
        if not error:
            #check if max/min angle not bigger than pi or -pi
            if MaxAngles[i*2+1] > np.around(np.pi,digitPi) or MinAngles[i*2+1] < -np.around(np.pi,digitPi):
                flagCountRel[actuatedHinges[i]-1,6]  += 1
                exflRel[(i+1)*stepsHinge-1] = 3
                error = True
            #check for max/min strech not bigger than a tolerance
            if (MaxStrRel[(i+1)*stepsHinge-1] > tolStretch or MinStrRel[(i+1)*stepsHinge-1]  < -tolStretch) and not error:
                flagCountRel[actuatedHinges[i]-1,7]  += 1
                exflRel[(i+1)*stepsHinge-1] = 4
                error = True
        if error:
            notConvHinges = np.append(notConvHinges,np.array([[i, exflFol[(i+1)*stepsHinge-1], exflRel[(i+1)*stepsHinge-1]]]), axis = 0)
            hingesMask[i] = np.NaN
#            print(hingeName[i])
            
    notConverged = sum(np.isnan(hingesMask))
    converged = totHingeNum - notConverged
    ############################################ normalize the flag counts
    allFlags = np.sum(np.add(flagCountFol,flagCountRel), axis = 0)
    allFlags = allFlags/totHingeNum
    for i in np.arange(totalflags):
        flagCountFol[:,i] = flagCountFol[:,i]/hingeCount
        flagCountRel[:,i] = flagCountRel[:,i]/hingeCount
        
    ############################################ ordering the last energies on the release according to the number of actuated hinges
#    lasteHingeRel = eHingeRel[stepsHinge-1::stepsHinge]
#    lasteEdgeRel = eEdgeRel[stepsHinge-1::stepsHinge]
##    for hinge in np.arange(totHingeNum):
##        if np.isnan(hingesMask[hinge]):
##            lasteHingeRel[hinge] = np.nan
##            lasteEdgeRel[hinge] = np.nan
#    lasteHingeRel = lasteHingeRel[orderedHinges]
#    lasteEdgeRel = lasteEdgeRel[orderedHinges]
#    energies = np.append(lasteHingeRel, lasteEdgeRel)
#    ############################################ going through all hinges in order and see all the neighbours in energy
#    print('Stable State hinge selection:')
#    differentEnergies = np.empty((0,2), dtype = int)
#    differentEnergiesName = np.empty(0)
#    differentEnergiesEnergy = np.empty((0,2))
#    for hinge in np.arange(totHingeNum):
#        if ~np.isnan(hingesMask[orderedHinges[hinge]]):
#            hingeeHinge = lasteHingeRel[hinge]
#            hingeeEdge = lasteEdgeRel[hinge]
#            sameEnergy = np.where(np.logical_and(np.logical_and(energies[:totHingeNum]>=hingeeHinge-tolHinge,
#                                                                energies[:totHingeNum]<=hingeeHinge+tolHinge),
#                                                 np.logical_and(energies[totHingeNum:]>=hingeeEdge-tolEdge, 
#                                                                energies[totHingeNum:]<=hingeeEdge+tolEdge)))[0]
#            ############################################ see if at least one of the neighbours converge (including the selected one)
#            i = 0
#            if len(notConvHinges) > 0:
#                #goes through all hinges in the same neighbourhood and see if they dont converge. 
#                #If one converges, it stops. If non converges gives i =-1
#                while len(np.where(notConvHinges[:,0] == orderedHinges[sameEnergy[i]])[0]) != 0:
#                    i = i + 1
#                    if i >= len(sameEnergy):
#                        i = -1
#                        break
#            ############################################ save the new state or not according to the convergance
#            if i != -1:
#                findit = np.where(differentEnergies[:,0] == orderedHinges[sameEnergy[i]])[0]
#                if len(findit) == 0:
#                    stablehinge = orderedHinges[hinge]
#                    differentEnergies = np.append(differentEnergies,np.array([[stablehinge, len(sameEnergy)]]), axis = 0)
#                    differentEnergiesName = np.append(differentEnergiesName, np.array(hingeName[stablehinge]))
#                    differentEnergiesEnergy = np.append(differentEnergiesEnergy, np.array([[lasteHingeRel[hinge],lasteEdgeRel[hinge]]]), axis = 0)
#                    print(hingeName[stablehinge])
#                elif findit >= 1:
#                    print(hingeName[orderedHinges[hinge]], findit) 
#            else:
#    #            nonConvStates += len(sameEnergy)
#                print('State from non convergent hinges')
#    #differentEnergies = np.array(differentEnergies)
    ######################### Problems: double counting of stable states that are close to each other due to the ill definition of neighbourhoods.
    #########################           Additionally counting states that didn't converge.
#    print('Total converged percentage:', converged/totHingeNum)
#    print('Not converged: ', notConverged, '/', totHingeNum)
    
    ###################### Analysis of angles to find stable states
    convHinges = np.empty(0, dtype = int)
    finalAngles = np.empty((0,np.size(dataAngles,1)))
    dataAngles = np.around(dataAngles/tolAngleSS)*tolAngleSS ## Here you conisder the tolerance for angles to recognize stable states
    for hinge in orderedHinges:
        if ~np.isnan(hingesMask[hinge]):
            sortAngleIndex = np.lexsort((dataAngles[2*hinge+1,:],dataAngles[2*hinge,:]))
            finalAngles = np.append(finalAngles, [dataAngles[2*hinge+1,sortAngleIndex]], axis = 0)
            convHinges = np.append(convHinges, hinge)
    
    differentAngles, index, counts = np.unique(finalAngles, axis = 0, return_index = True, return_counts = True)
    differentEnergies = np.column_stack((convHinges[index], counts))
    differentEnergiesName = hingeName[convHinges[index]]
    differentEnergiesEnergy = np.column_stack((eHingeRel[convHinges[index]*stepsHinge+stepsHinge-1], 
                                                eEdgeRel[convHinges[index]*stepsHinge+stepsHinge-1]))
   
    ###################### Normalize the angles to be proportional to Pi and shift the angle sum to easier understanding

    SumIntAngFol = SumIntAngFol/np.pi+internalHinges
    SumIntAngRel = SumIntAngRel/np.pi+internalHinges
    SumExtAngFol = (SumExtAngFol/np.pi-(totalnumberHinges-internalHinges))*(-1)
    SumExtAngRel = (SumExtAngRel/np.pi-(totalnumberHinges-internalHinges))*(-1)
    
    
    #%%
    #######################################################################################################################
    ##################### Ploting the result
    #######################################################################################################################
    if plot:
        fig1 = plt.figure(0,figsize=(cm2inch(35), cm2inch(20)))
        ax1 = plt.subplot(111)#
        NiceGraph2D(ax1, r'Average $\Delta$L',  r'Average $\Delta\theta$ [rad]')#, [min(eEdgeRel[stepsHinge-1::stepsHinge]), np.NaN],[0.0014,  np.NaN] )
    #    ax3 = plt.subplot(122)
    #    NiceGraph2D(ax3, 'Edge Energy',  'Hinge Energy', [0.00025, 0.012],[0.01375,  0.029] )
        
        fig2 = plt.figure(1,figsize=(cm2inch(35), cm2inch(20)))
        ax2 = plt.subplot(111)
        NiceGraph2D(ax2, 'Average Radius', 'StDev of Radius')#, [min(eEdgeRel[stepsHinge-1::stepsHinge]), np.NaN],[max(eEdgeRel[stepsHinge-1::stepsHinge]), np.NaN], buffer = [0.00001, 0.0])
        
        fig3 = plt.figure(2,figsize=(cm2inch(35), cm2inch(20)))
        ax4 = plt.subplot(111, projection='3d')
        NiceGraph3D(ax4, 'CenterMass X', 'CenterMass Y', 'CenterMass Z')
        
        fig4 = plt.figure(3,figsize=(cm2inch(35), cm2inch(20)))
        ax5 = plt.subplot(121)
        ax10 = plt.subplot(122)
        NiceGraph2D(ax5, r'Average $\Delta\theta$ [rad]', 'Sum of internal angles [rad]')
        NiceGraph2D(ax10, r'Average $\Delta\theta$ [rad]', 'Sum of external angles [rad]')
        ax5.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%g $\pi$'))
        ax10.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%g $\pi$'))
#        ax5.yaxis.set_major_locator(matl.ticker.MultipleLocator(base=1.0))
#        NiceGraph2D(ax5, 'Hinge-Set Number', 'Internal Hinge Energy', [np.nan, min(eHinIntRel)], [np.nan, max(eHinIntRel)], buffer = [0, 0.0004])
        
        fig5 = plt.figure(4,figsize=(cm2inch(35), cm2inch(20)))
        ax6 = plt.subplot(111)
#        ax6 = plt.subplot(121)  
#        ax7 = plt.subplot(122)   
        NiceGraph2D(ax6, '# of actuated Hinges', 'percentage # of flags', [0.5, 0], [len(hingeCount)+0.5, 1+0.01], [np.arange(len(hingeCount))+1, np.arange(0,1.1,0.1)])       
#        NiceGraph2D(ax7, '# of actuated Hinges', 'percentage # of flags', [0.5, 0], [len(hingeCount)+0.5, 1+0.01], [np.arange(len(hingeCount))+1, np.nan])  
        
        fig6 = plt.figure(5,figsize=(cm2inch(35), cm2inch(20)))
        ax8 = plt.subplot(111)
#        ax8 = plt.subplot(121)  
#        ax9 = plt.subplot(122) 
        NiceGraph2D(ax8, r'Average $\Delta$L', r'Max strecht $\Delta$L/L')
#        NiceGraph2D(ax9, 'Hinge-Set Number', 'Min final streching')            
                    
        width = 0.5
        colors = cm.Set2(np.linspace(0, 1, totalflags))
        for i, c in zip(reversed(np.arange(totalflags)), reversed(colors)):
            if np.sum(flagCountFol[:,i]+flagCountRel[:,i]) == 0:  ###### block to plot non-present flags
                continue
            ax6.bar(np.arange(len(hingeCount))+1, flagCountFol[:,i]+flagCountRel[:,i], width, color=c, label = i-3, bottom = np.sum(flagCountFol[:,:i],1)+np.sum(flagCountRel[:,:i],1))
#            ax6.bar(np.arange(len(hingeCount))+1, flagCountFol[:,i], width, color=c, bottom = np.sum(flagCountFol[:,:i],1), label = i-3)
#            ax7.bar(np.arange(len(hingeCount))+1, flagCountRel[:,i], width, color=c, bottom = np.sum(flagCountRel[:,:i],1))
        
        
        ax6.legend(loc = 2, fontsize =15, framealpha = 0.5, edgecolor = 'inherit', fancybox = False)
        s = 'Convergence: %.2f\nNon Convergence states: %d/%d' %(converged/totHingeNum, notConverged,totHingeNum)
        ax6.set_title(s)
        
        maxststs = np.max(differentEnergies[:,1])
        cmap1, norm1 = from_levels_and_colors(np.arange(np.nanmax(differentEnergies[:,1])+1)+1,
                                              cm.rainbow(np.linspace(0, 1, np.nanmax(differentEnergies[:,1]))))
        cmap1.set_over('r')
    
    
        for hinge in np.arange(totHingeNum):
            if ~np.isnan(hingesMask[hinge]):
                col = '#36648B'

            else:
                col = '#FE9128'
            findit = np.where(differentEnergies[:,0] == hinge)[0]
            if len(findit) != 0:
    #            ax3.scatter(eEdgeRel[stepsHinge*hinge+stepsHinge-1], eHingeRel[stepsHinge*hinge+stepsHinge-1], c = c)
                ax2.scatter(RadRel[stepsHinge*hinge+stepsHinge-1], StdRel[stepsHinge*hinge+stepsHinge-1], c = col, label = hingeName[hinge])
                ax4.scatter(CMxRel[stepsHinge*hinge+stepsHinge-1], CMyRel[stepsHinge*hinge+stepsHinge-1], CMzRel[stepsHinge*hinge+stepsHinge-1], c = col)
            ax5.scatter(eHingeRel[stepsHinge*hinge+stepsHinge-1], SumIntAngRel[stepsHinge*hinge+stepsHinge-1], c = col)
            ax10.scatter(eHingeRel[stepsHinge*hinge+stepsHinge-1], SumExtAngRel[stepsHinge*hinge+stepsHinge-1], c = col)
            ax8.scatter(eEdgeRel[stepsHinge*hinge+stepsHinge-1], abs(max(MaxStrRel[stepsHinge*hinge+stepsHinge-1],MinStrRel[stepsHinge*hinge+stepsHinge-1], key=abs)), c = col)
#                ax8.scatter(hingeNum[stepsHinge*hinge+stepsHinge-1], MaxStrRel[stepsHinge*hinge+stepsHinge-1], c = col)
#                ax9.scatter(hingeNum[stepsHinge*hinge+stepsHinge-1], abs(MinStrRel[stepsHinge*hinge+stepsHinge-1]), c = col)
#            ax1.plot(eEdgeRel[stepsHinge*hinge:stepsHinge*(hinge+1)],  eHingeRel[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = col)
#            ax1.scatter(eEdgeRel[stepsHinge*hinge+stepsHinge-1],  eHingeRel[stepsHinge*hinge+stepsHinge-1], c = col)                
    #            ax3.plot(eEdgeRel[stepsHinge*hinge:stepsHinge*(hinge+1)],  eHingeRel[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = col)
    #            ax2.plot(RadRel[stepsHinge*hinge:stepsHinge*(hinge+1)],  StdRel[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = col)
            if len(findit) != 0:# and differentEnergies[findit[0],1] > maxststs:
                ax1.annotate(hingeName[hinge], xy=(eEdgeRel[stepsHinge*hinge+stepsHinge-1], eHingeRel[stepsHinge*hinge+stepsHinge-1]), 
                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
                ax5.annotate(hingeName[hinge], xy=(eHingeRel[stepsHinge*hinge+stepsHinge-1], SumIntAngRel[stepsHinge*hinge+stepsHinge-1]), 
                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
                ax10.annotate(hingeName[hinge], xy=(eHingeRel[stepsHinge*hinge+stepsHinge-1], SumExtAngRel[stepsHinge*hinge+stepsHinge-1]), 
                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
                ax8.annotate(hingeName[hinge], xy=(eEdgeRel[stepsHinge*hinge+stepsHinge-1], abs(max(MaxStrRel[stepsHinge*hinge+stepsHinge-1],MinStrRel[stepsHinge*hinge+stepsHinge-1], key=abs))), 
                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
#                print(hingeName[differentEnergies[findit[0],0]], differentEnergies[findit[0],1])    
                
        cs1 = ax1.scatter(differentEnergiesEnergy[:,1], differentEnergiesEnergy[:,0], c = differentEnergies[:,1],
                          label = differentEnergiesName, cmap = cmap1, vmax = maxststs)
        
        cbar = plt.colorbar(cs1, format="%d", ax = ax1, fraction=0.05, pad=0.01, extend = 'max')
        cbar.set_ticks(np.linspace(1, maxststs, 13))
        cbar.set_label('Number of appearance', fontsize = 15, color = '0.2')
        cbar.ax.tick_params(axis='y',colors='0.2')
        cbar.ax.tick_params(axis='x',colors='0.2')
        cbar.outline.set_edgecolor('0.2')
                
        
        
    
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
        fig1.savefig(folder_name + file_name1[:-4]+normalized+'.png', transparent = False)
        #fig2.savefig(folder_name + file_name3[:-4]+normalized+'.png', transparent = True)
        #fig3.savefig(folder_name + '/CenterOfMass'+normalized+'.png', transparent = True)
        fig4.savefig(folder_name + '/InternalAnglesHingeEnergy'+normalized+'.png', transparent = True)
        fig5.savefig(folder_name + '/Flags.png', transparent = True)
        fig6.savefig(folder_name + '/MaxStretchEdgeEnergy'+normalized+'.png', transparent = True)



    return allFlags, len(differentEnergies[:,0]), differentEnergiesName, differentEnergiesEnergy, differentAngles, SumIntAngRel[stepsHinge*index+stepsHinge-1], SumExtAngRel[stepsHinge*index+stepsHinge-1]

#hinges = np.zeros(30)
#
#with open(folder_name + file_name1,'r',newline='') as ffit: 
#    readers = csv.reader(ffit,delimiter=',');
#    for row in readers:
#        if row:
#            hinges[np.size(row)] +=1
#%%
if __name__ == "__main__":
    folder_name = "Results/cube/sqp/energy/15-Dec-2017_100AngleCnstr\kh0.001_kta0.562_ke0.316_kf1.000"
    ReadandAnalizeFile(folder_name, khinge = 0.001, kedge = 10)