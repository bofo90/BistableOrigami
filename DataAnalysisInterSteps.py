import matplotlib as matl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import from_levels_and_colors
import numpy as np
import configparser
import os.path
import scipy.cluster.hierarchy as hierarch
from scipy.spatial.distance import pdist
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

def ReadandAnalizeFile(folder_name, plot = False, normalize = False):
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
    eEdge = dataEnergy[1,:]
    eDiag = dataEnergy[2,:]
    eFace = dataEnergy[3,:]
    eHinge = dataEnergy[4,:]
    eTAngle = dataEnergy[5,:]
    eHinInt = dataEnergy[6,:]
    exfl = dataEnergy[7,:].astype(int)
    
    eAllEdge = eEdge + eFace
        
    hingeName = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [1], dtype=bytes).astype(str)    
        
    dataPosStad = np.loadtxt(folder_name+file_name3,skiprows=1, delimiter = ',', unpack = True)
    CMx = dataPosStad[1,:]
    CMy = dataPosStad[2,:]
    CMz = dataPosStad[3,:]
    Rad = dataPosStad[4,:]
    Std = dataPosStad[5,:]
    MaxStr = dataPosStad[6,:]
    MinStr = dataPosStad[7,:]
    SumIntAng = dataPosStad[8,:]
    SumExtAng = dataPosStad[9,:]
    
    dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
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
        totalnumberDiag = int(metadata.get('extUnitCell', 'Diag'))
        khinge = float(metadata.get('options','kHinge'))
        kedge = float(metadata.get('options','kEdge'))
        kdiag = float(metadata.get('options','kDiag'))
        kface = float(metadata.get('options','kFace'))
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 
        #### If there is no metadata file the order of the dataPosStad is not the same as here shown           
        
    tolStretch = 0.3 #max precentage of allowed stretch
    digitPi = 4 # digits of pi for rounding to see if angles go beyond it and "not converge"
        
    stepsHinge = int(len(hingeNum)/len(hingeName))
    totHingeNum = len(hingeName)
    totalflags = 11
    #%%
    #######################################################################################################################
    ##################### Analize the data
    #######################################################################################################################


    ############################################ Modify data to have a normalized energy (or just difference in angle/length)
    if normalize:
        eHinge = np.sqrt(eHinge*2/khinge/totalnumberHinges)
        eHinInt = np.sqrt(eHinInt*2/khinge/internalHinges)
        eEdge = eEdge*2/kedge
        eDiag = eDiag*2/kdiag
        eAllEdge = np.sqrt((eEdge + eDiag)/(totalnumberEdges+totalnumberDiag))
        eEdge = np.sqrt(eEdge/totalnumberEdges)
        eDiag = np.sqrt(eDiag/totalnumberDiag)
        normalized = 'norm'
    
    #%%
    ############################################ get the number of actuated hinges for each hinge-set
    actuatedHinges = np.zeros(totHingeNum, dtype = int)
    hingeCount = np.zeros(internalHinges)   
    allHinges = np.arange(totHingeNum)
    for hinge in allHinges:
        actuatedHinges[hinge] = len(hingeName[hinge].split())
        hingeCount[actuatedHinges[hinge]-1] += 1
        
    #%%
    ############################################ order the hinge-set according to the number of actuated hinges
    orderedHinges = np.argsort(actuatedHinges)
    
    ############################################ count the different flags and Mask the hinge-sets that didnt converge
    exfl = exfl.reshape(totHingeNum,stepsHinge)
    flagmask = np.logical_and(exfl !=1,exfl !=5)
    
    
    
    flagCountFol = np.zeros((internalHinges, totalflags))
    flagCountRel = np.zeros((internalHinges, totalflags))
    notConvHinges = np.empty((0,3), dtype = int)
    for i in allHinges:
        if flagmask[i,0]:
            flagCountFol[actuatedHinges[i]-1,exfl[i,0]]  += 1
        if flagmask[i,1]:# and exflRel[(i+1)*stepsHinge-1] != -2:
            flagCountRel[actuatedHinges[i]-1,exfl[i,1]]  += 1
        if flagmask[i,:].any():
            #check if max/min angle not bigger than pi or -pi
            if MaxAngles[i*3+2] > np.around(np.pi,digitPi) or MinAngles[i*3+2] < -np.around(np.pi,digitPi):
                flagCountRel[actuatedHinges[i]-1,6]  += 1
                exfl[i,1] = 6
                flagmask[i,1] = True
            #check for max/min strech not bigger than a tolerance
            if (MaxStr[(i+1)*stepsHinge-1] > tolStretch or MinStr[(i+1)*stepsHinge-1]  < -tolStretch) and not flagmask[i,1]:
                flagCountRel[actuatedHinges[i]-1,7]  += 1
                exfl[i,1] = 7
                flagmask[i,1] = True
        if flagmask[i,:].any():
            notConvHinges = np.append(notConvHinges,np.array([[i, exfl[i,0],exfl[i,1]]]), axis = 0)
#            hingesMask[i] = np.NaN
#            print(hingeName[i])
            
    notConverged = sum(flagmask.any(axis = 1))
    converged = totHingeNum - notConverged
    hingesMask = np.logical_not(flagmask.any(axis = 1))
    convHinges = allHinges[hingesMask]
    ############################################ normalize the flag counts
    allFlags = np.sum(np.add(flagCountFol,flagCountRel), axis = 0)
    allFlags = allFlags/totHingeNum
    for i in np.arange(totalflags):
        flagCountFol[:,i] = flagCountFol[:,i]/hingeCount
        flagCountRel[:,i] = flagCountRel[:,i]/hingeCount
        
    #%%
    ###################### Normalize the angles to be proportional to Pi and shift the angle sum to easier understanding
    SumIntAng = SumIntAng/np.pi+internalHinges
    SumExtAng = (SumExtAng/np.pi-(totalnumberHinges-internalHinges))*(-1)

    #%%
    ###################### Analysis of angles to find stable states
    
    finalAngles = np.empty((0,np.size(dataAngles,1)))
    for hinge in np.arange(totHingeNum):
        sortAllAngIndex = np.lexsort((dataAngles[(stepsHinge+1)*(hinge+1)-1,:],dataAngles[(stepsHinge+1)*hinge,:]))
        finalAngles = np.append(finalAngles, [dataAngles[(stepsHinge+1)*(hinge+1)-1,sortAllAngIndex]], axis = 0)
        
    Z = hierarch.linkage(finalAngles[hingesMask], 'centroid')
    inverse = hierarch.fcluster(Z, 1, criterion='distance')
    c = hierarch.cophenet(Z, pdist(finalAngles[hingesMask]))
    print('this is the cophenet of the hierarchical linkage', c[0])
    
    ###################### Plot the cluster and see how are the results related
#    plt.figure(figsize=(25, 10))
#    plt.title('Hierarchical Clustering Dendrogram')
#    plt.xlabel('sample index')
#    plt.ylabel('distance')
#    hierarch.dendrogram(
#        Z,
#        truncate_mode='lastp',  # show only the last p merged clusters
#        p=50,  # show only the last p merged clusters
#        leaf_rotation=90.,  # rotates the x axis labels
#        leaf_font_size=8.,  # font size for the x axis labels
#        show_contracted=True,
#    )
#    plt.show()

    
    SS, SSposmasked, SScounts = np.unique(inverse, return_index = True, return_counts=True)
    SSpos = convHinges[SSposmasked]
    differentEnergies = np.column_stack((SSpos, SScounts))
    differentEnergiesName = hingeName[SSpos]
    differentEnergiesEnergy = np.column_stack((eHinge[SSpos*stepsHinge+stepsHinge-1], 
                                                eAllEdge[SSpos*stepsHinge+stepsHinge-1]))
    


    if np.size(SS) == 1:
        print('Error: No additional stable states found.\n')
        if plot:
            fig5 = plt.figure(4,figsize=(cm2inch(35), cm2inch(20)))
            ax6 = plt.subplot(111)
            NiceGraph2D(ax6, '# of actuated Hinges', 'percentage # of flags', [0.5, 0], [len(hingeCount)+0.5, 1+0.01], [np.arange(len(hingeCount))+1, np.arange(0,1.1,0.1)])       
            width = 0.5
            colors = cm.tab10(np.linspace(0, 1, totalflags))
            for i, c in zip(reversed(np.arange(totalflags)), reversed(colors)):
                if np.sum(flagCountFol[:,i]+flagCountRel[:,i]) == 0:  ###### block to plot non-present flags
                    continue
                ax6.bar(np.arange(len(hingeCount))+1, flagCountFol[:,i]+flagCountRel[:,i], width, color=c, label = i-3, bottom = np.sum(flagCountFol[:,:i],1)+np.sum(flagCountRel[:,:i],1))
            ax6.legend(loc = 2, fontsize =15, framealpha = 0.5, edgecolor = 'inherit', fancybox = False)
            s = 'Convergence: %.2f\nNon Convergence states: %d/%d' %(converged/totHingeNum, notConverged,totHingeNum)
            ax6.set_title(s)
        return allFlags, 0, [], [], [], [], []

    
    #%%
    #######################################################################################################################
    ##################### Ploting the result
    #######################################################################################################################
    if plot:
        fig1 = plt.figure(0,figsize=(cm2inch(35), cm2inch(20)))
        ax1 = plt.subplot(111)#
        NiceGraph2D(ax1, r'Average $\Delta$L',  r'Average $\Delta\theta$ [rad]')#, [min(eAllEdgeRel[stepsHinge-1::stepsHinge]), np.NaN],[0.0014,  np.NaN] )
    #    ax3 = plt.subplot(122)
    #    NiceGraph2D(ax3, 'Edge Energy',  'Hinge Energy', [0.00025, 0.012],[0.01375,  0.029] )
        
        fig2 = plt.figure(1,figsize=(cm2inch(35), cm2inch(20)))
        ax2 = plt.subplot(111)
        NiceGraph2D(ax2, 'Average Radius', 'StDev of Radius')#, [min(eAllEdgeRel[stepsHinge-1::stepsHinge]), np.NaN],[max(eAllEdgeRel[stepsHinge-1::stepsHinge]), np.NaN], buffer = [0.00001, 0.0])
        
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
        colors = cm.tab10(np.linspace(0, 1, totalflags))
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
    
    
        for hinge in allHinges:
            if hingesMask[hinge]:
                col = '#36648B'
            else:
                col = '#FE9128'
            findit = np.where(differentEnergies[:,0] == hinge)[0]
            i = hinge*stepsHinge+stepsHinge-1
            if len(findit) != 0:
#                ax3.scatter(eAllEdge[i], eHinge[i], c = c)
                ax2.scatter(Rad[i], Std[i], c = col, label = hingeName[hinge])
                ax4.scatter(CMx[i], CMy[i], CMz[i], c = col)
                ax5.scatter(eHinge[i], SumIntAng[i], c = col)
                ax10.scatter(eHinge[i], SumExtAng[i], c = col)
                ax8.scatter(eAllEdge[i], abs(max(MaxStr[i],MinStr[i], key=abs)), c = col)
#                ax8.scatter(hingeNum[i], MaxStr[i], c = col)
#                ax9.scatter(hingeNum[i], abs(MinStrRel[i]), c = col)
#            ax1.plot(eAllEdge[stepsHinge*hinge:stepsHinge*(hinge+1)],  eHinge[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = col)
#            ax1.scatter(eAllEdge[i],  eHinge[i], c = col)                
    #            ax3.plot(eAllEdge[stepsHinge*hinge:stepsHinge*(hinge+1)],  eHinge[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = col)
    #            ax2.plot(Rad[stepsHinge*hinge:stepsHinge*(hinge+1)],  Std[stepsHinge*hinge:stepsHinge*(hinge+1)], '--',c = col)
            if len(findit) != 0:# and differentEnergies[findit[0],1] > maxststs:
#                ax1.annotate(hingeName[hinge], xy=(eAllEdge[i], eHinge[i]), 
#                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
#                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
#                ax5.annotate(hingeName[hinge], xy=(eHinge[i], SumIntAng[i]), 
#                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
#                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
#                ax10.annotate(hingeName[hinge], xy=(eHinge[i], SumExtAng[i]), 
#                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
#                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
#                ax8.annotate(hingeName[hinge], xy=(eAllEdge[i], abs(max(MaxStr[i],MinStr[i], key=abs))), 
#                              xytext=(10, 10), textcoords='offset points', ha='right', va='bottom',
#                              arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
                print(hingeName[differentEnergies[findit[0],0]], differentEnergies[findit[0],1])    
                
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



    return allFlags, len(differentEnergies[:,0]), differentEnergiesName, differentEnergiesEnergy, differentAngles, SumIntAngRel[stepsHinge*convHinges[index]+stepsHinge-1], SumExtAngRel[stepsHinge*convHinges[index]+stepsHinge-1]

#hinges = np.zeros(30)
#
#with open(folder_name + file_name1,'r',newline='') as ffit: 
#    readers = csv.reader(ffit,delimiter=',');
#    for row in readers:
#        if row:
#            hinges[np.size(row)] +=1
#%%
if __name__ == "__main__":
    folder_name = "Results/tetrahedron/active-set/energy/13-Jun-2018_temp\kh0.001_kta100.000_ke3.000_kf100.000"
    ReadandAnalizeFile(folder_name)