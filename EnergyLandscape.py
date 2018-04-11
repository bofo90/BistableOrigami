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
#import seaborn as sns

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
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
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

folder_name = "Results/truncated tetrahedron/active-set/energy/15-Feb-2018_EnergyAllAngles_3_24\kh0.001_kta100.000_ke3.000_kf100.000"
inverted = False
maxEnergy = 0.75
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
exfl = dataEnergy[6,:].astype(int)

hingeName = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [1], dtype=bytes).astype(str)   
closingAngl1 = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [2])/np.pi
closingAngl2 = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [3])/np.pi
angleNum = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', usecols = [4,5])

dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
dataAngles = np.delete(dataAngles, 0, 1)
MaxAngles = np.max(dataAngles, axis = 1)
MinAngles = np.min(dataAngles, axis = 1)

TotSimul = np.size(hingeName)
IterPerSimul = np.int(np.size(dataAngles,0)/TotSimul)
IterPerSimulEnergy = np.int(np.size(dataEnergy,1)/TotSimul)
NumberAngles = np.size(dataAngles,1)

FoldAng = np.array(hingeName[0][1:-1].split(), int)
FoldAng = np.sort(FoldAng)

metadataFile = '/metadata.txt'

if os.path.isfile(folder_name+metadataFile):
    metadata = configparser.RawConfigParser()
    metadata.read(folder_name+metadataFile)

    totalnumberHinges = int(metadata.get('extUnitCell', 'Hinges'))
    totalnumberEdges = int(metadata.get('extUnitCell', 'Edges'))
    totalnumberDiag = int(metadata.get('extUnitCell', 'Diag'))
    khinge = float(metadata.get('options','kHinge'))
    kedge = float(metadata.get('options','kEdge'))
    kdiag = float(metadata.get('options','kDiag'))
    kface = float(metadata.get('options','kFace'))
else:
    raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

eEdge = eEdge/kedge/totalnumberEdges
eDiag = eDiag/kdiag/totalnumberDiag
eHinge = eHinge/khinge/totalnumberHinges

eTotal= eEdge + eDiag + eHinge

#%%
###Check if a folding pattern didnt converge
exfl = exfl.reshape(TotSimul,IterPerSimulEnergy)
flagmask = np.logical_and(exfl !=1,exfl !=5)
flagmask = flagmask.any(axis= 1)
if flagmask.any():
    print('Error: There was at least one non convergent fold pattern.\n')

##Do the energy landscape plot
divitheta1 = len(np.unique(closingAngl1))
divitheta2 = len(np.unique(closingAngl2))

tickstheta1 = 5
tickstheta2 = 5

sortAngl = np.lexsort((closingAngl2, closingAngl1))
closingAngl1 = -closingAngl1[sortAngl[::-1]]
closingAngl2 = -closingAngl2[sortAngl[::-1]]

theta1 = closingAngl1.reshape((divitheta1,divitheta2))
theta2 = closingAngl2.reshape((divitheta1,divitheta2))

totEnergysort = eTotal[IterPerSimulEnergy-2::IterPerSimulEnergy]
#totEnergysort = eDiag[IterPerSimulEnergy-2::IterPerSimulEnergy]+eEdge[IterPerSimulEnergy-2::IterPerSimulEnergy]
#totEnergysort = eEdge[IterPerSimulEnergy-2::IterPerSimulEnergy]
#totEnergysort = eDiag[IterPerSimulEnergy-2::IterPerSimulEnergy]
#totEnergysort = eHinge[IterPerSimulEnergy-2::IterPerSimulEnergy]
#totEnergysort = eTAngle[IterPerSimulEnergy-2::IterPerSimulEnergy]

totEnergysort = np.ma.masked_array(totEnergysort, mask=flagmask)
totEnergysort = totEnergysort[sortAngl[::-1]]
if inverted:
    totEnergyMat = totEnergysort.reshape((divitheta1,divitheta2))
else:
    totEnergyMat = totEnergysort.reshape((divitheta1,divitheta2)).T

sep1 = (closingAngl1[-1]-closingAngl1[0])/(divitheta1-1)/2
sep2 = (closingAngl2[-1]-closingAngl2[0])/(divitheta2-1)/2

fig1 = plt.figure(0,figsize=(cm2inch(24.1), cm2inch(20)))
ax1 = plt.subplot(111)
if inverted:
    NiceGraph2D(ax1, 'TargAngl Hinge %d [rad]' %FoldAng[0], 'TargAngl Hinge %d [rad]'%FoldAng[1], mincoord = [0,0], maxcoord = [1,1],divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
                
#                mincoord = [closingAngl2[0], closingAngl1[0]], 
#                maxcoord = [closingAngl2[-1], closingAngl1[-1]],  divisions = [tickstheta2, tickstheta1], buffer = [sep2, sep1])
else:
    NiceGraph2D(ax1, 'TargAngl Hinge %d [rad]' %FoldAng[0], 'TargAngl Hinge %d [rad]'%FoldAng[1], mincoord = [0,0], maxcoord = [1,1],divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
                
#                mincoord = [closingAngl1[0], closingAngl2[0]], 
#                maxcoord = [closingAngl1[-1], closingAngl2[-1]],  divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])

if not inverted:
    cs1 = ax1.imshow(totEnergyMat, extent=[theta2[0,0]-sep2,theta2[0,-1]+sep2,theta1[0,0]-sep1,theta1[-1,0]+sep1], 
                     cmap = cm.nipy_spectral, aspect = 'auto',vmax = maxEnergy, origin = 'lower') #nipy_spectral
else:
    cs1 = ax1.imshow(totEnergyMat, extent=[theta1[0,0]-sep1,theta1[-1,0]+sep1,theta2[0,0]-sep2,theta2[0,-1]+sep2], 
                     cmap = cm.nipy_spectral, aspect = 'auto',vmax = maxEnergy, origin = 'lower')

SStheta1 = -dataAngles[3::IterPerSimul,FoldAng[0]-1]/np.pi
SStheta1 = SStheta1[sortAngl[::-1]]
SStheta2 = -dataAngles[3::IterPerSimul,FoldAng[1]-1]/np.pi
SStheta2 = SStheta2[sortAngl[::-1]]

ax1.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))
ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))

cbar = plt.colorbar(cs1,  ax = ax1, fraction=0.05, pad=0.01, extend = 'max')#format="%d", 
cbar.set_ticks(np.linspace(0, maxEnergy, 5))
cbar.set_label('Energy', fontsize = 15, color = '0.2')
cbar.ax.tick_params(axis='y',colors='0.2')
cbar.ax.tick_params(axis='x',colors='0.2')
cbar.outline.set_edgecolor('0.2')

######################################################################
#Analysis for stable states

finalAngles = np.empty((0,np.size(dataAngles,1)))
for hinge in sortAngl[::-1]:
    sortAllAngIndex = np.lexsort((dataAngles[IterPerSimul*(hinge+1)-1,:],dataAngles[IterPerSimul*hinge,:]))
    finalAngles = np.append(finalAngles, [dataAngles[IterPerSimul*(hinge+1)-1,sortAllAngIndex]], axis = 0)


Z = hierarch.linkage(finalAngles, 'centroid')
inverse = hierarch.fcluster(Z, 1, criterion='distance')
c = hierarch.cophenet(Z, pdist(finalAngles))
print('this is the cophenet of the hierarchical linkage', c[0])

SS, SSpos = np.unique(inverse, return_index = True)
angleNum = angleNum[sortAngl[::-1],:]
print(angleNum[SSpos,:])

inverse = np.ma.masked_array(inverse, mask=flagmask[sortAngl[::-1]])

if inverted:
    stableStateMat = inverse.reshape((divitheta1,divitheta2))
else:
    stableStateMat = inverse.reshape((divitheta1,divitheta2)).T

fig2 = plt.figure(1,figsize=(cm2inch(24.1), cm2inch(20)))
ax2 = plt.subplot(111)

if not inverted:
    NiceGraph2D(ax2, 'TargAngl Hinge %d [rad]' %FoldAng[0], 'TargAngl Hinge %d [rad]'%FoldAng[1], 
                mincoord = [0,0], maxcoord = [1,1],divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
#                mincoord = [closingAngl2[0], closingAngl1[0]], 
#                maxcoord = [closingAngl2[-1], closingAngl1[-1]],  divisions = [tickstheta2,tickstheta1], buffer = [sep2, sep1])
else:
    NiceGraph2D(ax2, 'TargAngl Hinge %d [rad]' %FoldAng[0], 'TargAngl Hinge %d [rad]'%FoldAng[1], 
                mincoord = [0,0], maxcoord = [1,1],divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
#                mincoord = [closingAngl1[0], closingAngl2[0]], 
#                maxcoord = [closingAngl1[-1], closingAngl2[-1]],  divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])

cmap2, norm2 = from_levels_and_colors(np.linspace(0,np.max(inverse),np.max(inverse)+1),
                                      cm.Set3(np.linspace(0, 1, np.max(inverse)))) #gist_rainbow #Set3

if not inverted:
    cs2 = ax2.imshow(stableStateMat, extent=[theta2[0,0]-sep2,theta2[0,-1]+sep2,theta1[0,0]-sep1,theta1[-1,0]+sep1], 
                     cmap = cmap2, aspect = 'auto', origin = 'lower')
else:
    cs2 = ax2.imshow(stableStateMat, extent=[theta1[0,0]-sep1,theta1[-1,0]+sep1,theta2[0,0]-sep2,theta2[0,-1]+sep2], 
                     cmap = cmap2, aspect = 'auto', origin = 'lower')

ax2.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))
ax2.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))

cbar2 = plt.colorbar(cs2,  ax = ax2, fraction=0.05, pad=0.01)#, extend = 'max'
cbar2.set_ticks(np.linspace(0, np.max(inverse), np.max(inverse)+1))
cbar2.set_label('Stable State', fontsize = 15, color = '0.2')
cbar2.ax.tick_params(axis='y',colors='0.2')
cbar2.ax.tick_params(axis='x',colors='0.2')
cbar2.outline.set_edgecolor('0.2')

#############################################################################
#Plot of real angles

realtheta1 = -dataAngles[2::IterPerSimul,FoldAng[0]-1]/np.pi
realtheta1 = realtheta1[sortAngl[::-1]]
realtheta2 = -dataAngles[2::IterPerSimul,FoldAng[1]-1]/np.pi
realtheta2 = realtheta2[sortAngl[::-1]]

fig3 = plt.figure(2,figsize=(cm2inch(24.1), cm2inch(20)))
ax3 = plt.subplot(111)
if inverted:
    NiceGraph2D(ax3, 'Hinge %d [rad]' %FoldAng[0], 'Hinge Hinge %d [rad]'%FoldAng[1], 
                mincoord = [0,0], maxcoord = [1,1],divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
#                mincoord = [min(realtheta1), min(realtheta2)], 
#                maxcoord = [max(realtheta1), max(realtheta2)],  divisions = [tickstheta2, tickstheta1], buffer = [sep2, sep1])
else:
    NiceGraph2D(ax3, 'Hinge %d [rad]' %FoldAng[0], 'Hinge Hinge %d [rad]'%FoldAng[1],
                mincoord = [0,0], maxcoord = [1,1],divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])
#                mincoord = [min(realtheta1), min(realtheta2)], 
#                maxcoord = [max(realtheta1), max(realtheta2)],  divisions = [tickstheta1, tickstheta2], buffer = [sep1, sep2])

cmap3, norm3 = from_levels_and_colors(np.linspace(0,maxEnergy,1000), cm.rainbow(np.linspace(0, 1, 1000)), extend = 'max')

cs3 = ax3.scatter(realtheta1, realtheta2, c = totEnergysort, cmap = cmap3, vmax = maxEnergy, s = 96, marker = 's') #150 #360

ax3.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))
ax3.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))

cbar3 = plt.colorbar(cs3, ax = ax3, fraction=0.05, pad=0.01, extend = 'max')
cbar3.set_ticks(np.linspace(0, maxEnergy, 5))
cbar3.set_label('Energy', fontsize = 15, color = '0.2')
cbar3.ax.tick_params(axis='y',colors='0.2')
cbar3.ax.tick_params(axis='x',colors='0.2')
cbar3.outline.set_edgecolor('0.2')

#############################################################################
#Modifications of plots

#adding boundaries of stables states on energy landscape plot

if not inverted:
    for line in np.arange(divitheta1):
        for row in np.arange(divitheta2):
            #Vertical lines of separation
            if line != divitheta1-1:
                if not (np.ma.is_masked(stableStateMat[row,line]) or np.ma.is_masked(stableStateMat[row,line+1])):
                    if stableStateMat[row,line] != stableStateMat[row,line+1]:
    #                    ax1.plot([theta2[row,line]+sep2,theta2[row,line]+sep2],[theta1[row,line]-sep1,theta1[row,line]+sep1] ,c='k', linewidth = 1)
                        ax1.plot([theta2[0,0]+sep2+(2*sep2)*(line),theta2[0,0]+sep2+(2*sep2)*(line)],
                                  [theta1[0,0]-sep1+(2*sep1)*row,theta1[0,0]-sep1+(2*sep1)*(row+1)] ,c='k', linewidth = 1)
    
            #Horizontal lines of separation
            if row != divitheta2-1:
                if not (np.ma.is_masked(stableStateMat[row,line]) or np.ma.is_masked(stableStateMat[row+1,line])):
                    if stableStateMat[row+1,line] != stableStateMat[row,line]:
    #                    ax1.plot([theta2[row,line]-sep2,theta2[row,line]+sep2],[theta1[row,line]+sep1,theta1[row,line]+sep1] ,c='k', linewidth = 1.5)
                        ax1.plot([theta2[0,0]-sep2+(2*sep2)*(line),theta2[0,0]-sep2+(2*sep2)*(line+1)],
                                  [theta1[0,0]+sep1+(2*sep1)*row,theta1[0,0]+sep1+(2*sep1)*(row)] ,c='k', linewidth = 1)
else:
    for line in np.arange(divitheta2):
        for row in np.arange(divitheta1):
            #Vertical lines of separation
            if line != divitheta2-1:
                if not (np.ma.is_masked(stableStateMat[row,line]) or np.ma.is_masked(stableStateMat[row,line+1])):
                    if stableStateMat[row,line] != stableStateMat[row,line+1]:
    #                    ax1.plot([theta2[row,line]+sep2,theta2[row,line]+sep2],[theta1[row,line]-sep1,theta1[row,line]+sep1] ,c='k', linewidth = 1)
                        ax1.plot([theta1[0,0]+sep1+(2*sep2)*(line),theta1[0,0]+sep1+(2*sep1)*(line)],
                                  [theta2[0,0]-sep2+(2*sep2)*row,theta2[0,0]-sep2+(2*sep2)*(row+1)] ,c='k', linewidth = 1)
    
            #Horizontal lines of separation
            if row != divitheta1-1:
                if not (np.ma.is_masked(stableStateMat[row,line]) or np.ma.is_masked(stableStateMat[row+1,line])):
                    if stableStateMat[row+1,line] != stableStateMat[row,line]:
    #                    ax1.plot([theta2[row,line]-sep2,theta2[row,line]+sep2],[theta1[row,line]+sep1,theta1[row,line]+sep1] ,c='k', linewidth = 1.5)
                        ax1.plot([theta1[0,0]-sep1+(2*sep1)*(line),theta1[0,0]-sep1+(2*sep1)*(line+1)],
                                  [theta2[0,0]+sep2+(2*sep2)*row,theta2[0,0]+sep2+(2*sep2)*(row)] ,c='k', linewidth = 1)
          
#adding stars on final angles of stable states
ax1.scatter(SStheta1, SStheta2, c = inverse, cmap = cmap2, s = 200, marker = '*', edgecolor = 'k', lw = 0.2, zorder = 3)  
ax2.scatter(SStheta1, SStheta2, c = inverse, cmap = cmap2, s = 250, marker = '*', edgecolor = 'k', lw = 0.2)              
            

fig1.tight_layout()
fig1.show()
fig1.savefig(folder_name + '/EnergyLand.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyAllEdges.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyEdge.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyDiag.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyHinge.png', transparent = True)
#fig1.savefig(folder_name + '/EnergyTA.png', transparent = True)

fig2.tight_layout()
fig2.show()
fig2.savefig(folder_name + '/StableStates.png', transparent = True)

fig3.tight_layout()
fig3.show()
fig3.savefig(folder_name + '/RealAngles-Energy.png', transparent = True)