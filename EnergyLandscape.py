import matplotlib as matl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import from_levels_and_colors
import numpy as np
import configparser
import os.path


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

folder_name = "Results/truncated tetrahedron/active-set/energy/08-Feb-2018_Energylandscape_3to24_b\kh0.001_kta100.000_ke1.000_kf100.000"
inverted = False
tolAngleSS = 0.174 # equivalent to 10 degrees
maxEnergy = 0.44
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
eDiagFol = dataEnergy[3,:]
eDiagRel = dataEnergy[4,:]    
eFaceFol = dataEnergy[5,:]
eFaceRel = dataEnergy[6,:]
eHingeFol = dataEnergy[7,:]
eHingeRel = dataEnergy[8,:]
eTAngleFol = dataEnergy[9,:]
eTAngleRel = dataEnergy[10,:]
eHinIntFol = dataEnergy[11,:]
eHinIntRel = dataEnergy[12,:]
exflFol = dataEnergy[13,:].astype(int)
exflRel = dataEnergy[14,:].astype(int)

eTotalFol= eEdgeFol + eDiagFol +eFaceFol + eHingeFol
eTotalRel= eEdgeRel + eDiagRel +eFaceRel + eHingeRel

hingeName = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [1], dtype=bytes).astype(str)   
closingAngl1 =  np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [2])/np.pi
closingAngl2 =  np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', unpack = True, usecols = [3])/np.pi

dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
dataAngles = np.delete(dataAngles, 0, 1)
MaxAngles = np.max(dataAngles, axis = 1)
MinAngles = np.min(dataAngles, axis = 1)

TotSimul = np.size(closingAngl1)
IterPerSimul = np.int(np.size(dataAngles,0)/TotSimul)

#%%
###Check if a folding pattern didnt converge
if len(np.unique(exflFol))>2:
    print('Error: There was at least one non convergent fold pattern.\n')
if len(np.unique(exflRel))>2:
    print('Error: There was at least one non convergent release pattern.\n')


##Do the energy landscape plot
divitheta1 = len(np.unique(closingAngl1))
divitheta2 = len(np.unique(closingAngl2))

sortAngl = np.lexsort((closingAngl2, closingAngl1))
closingAngl1 = -closingAngl1[sortAngl[::-1]]
closingAngl2 = -closingAngl2[sortAngl[::-1]]

theta1 = closingAngl1.reshape((divitheta1,divitheta2))
theta2 = closingAngl2.reshape((divitheta1,divitheta2))

totEnergysort = eTotalFol[1::2]
#totEnergysort = eTAngleFol[1::2]
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
    NiceGraph2D(ax1, 'Hinge 3 [rad]', 'Hinge 24 [rad]',mincoord = [closingAngl2[0], closingAngl1[0]], 
                maxcoord = [closingAngl2[-1], closingAngl1[-1]],  divisions = [divitheta2, divitheta1], buffer = [sep2, sep1])
else:
    NiceGraph2D(ax1, 'Hinge 3 [rad]', 'Hinge 24 [rad]',mincoord = [closingAngl1[0], closingAngl2[0]], 
                maxcoord = [closingAngl1[-1], closingAngl2[-1]],  divisions = [divitheta1, divitheta2], buffer = [sep1, sep2])

if inverted:
    cs1 = ax1.imshow(totEnergyMat, extent=[theta2[0,0]-sep2,theta2[0,-1]+sep2,theta1[0,0]-sep1,theta1[-1,0]+sep1], 
                     cmap = cm.rainbow, aspect = 'auto',vmax = maxEnergy, origin = 'lower')
else:
    cs1 = ax1.imshow(totEnergyMat, extent=[theta1[0,0]-sep1,theta1[-1,0]+sep1,theta2[0,0]-sep2,theta2[0,-1]+sep2], 
                     cmap = cm.rainbow, aspect = 'auto',vmax = maxEnergy, origin = 'lower')

ax1.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))
ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))

cbar = plt.colorbar(cs1,  ax = ax1, fraction=0.05, pad=0.01, extend = 'max')#format="%d", 
cbar.set_ticks(np.linspace(0, maxEnergy, 5))
cbar.set_label('Energy', fontsize = 15, color = '0.2')
cbar.ax.tick_params(axis='y',colors='0.2')
cbar.ax.tick_params(axis='x',colors='0.2')
cbar.outline.set_edgecolor('0.2')

fig1.tight_layout()
fig1.show()
fig1.savefig(folder_name + '/EnergyLand.png', transparent = False)
#fig1.savefig(folder_name + '/EnergyTA.png', transparent = False)

######################################################################
#Analysis for stable states

finalAngles = np.empty((0,np.size(dataAngles,1)))
dataAnglesNorm = np.around(dataAngles/tolAngleSS)*tolAngleSS ## Here you conisder the tolerance for angles to recognize stable states
for hinge in sortAngl[::-1]:
    sortAllAngIndex = np.lexsort((dataAnglesNorm[IterPerSimul*(hinge+1)-1,:],dataAnglesNorm[IterPerSimul*hinge,:]))
    finalAngles = np.append(finalAngles, [dataAnglesNorm[IterPerSimul*(hinge+1)-1,sortAllAngIndex]], axis = 0)

differentAngles, index, inverse, counts = np.unique(finalAngles, axis = 0, return_index = True, return_inverse = True, return_counts = True)
differentEnergies = np.column_stack((sortAngl[index], counts))

if inverted:
    stableStateMat = inverse.reshape((divitheta1,divitheta2))
else:
    stableStateMat = inverse.reshape((divitheta1,divitheta2)).T

fig2 = plt.figure(1,figsize=(cm2inch(24.1), cm2inch(20)))
ax2 = plt.subplot(111)

if inverted:
    NiceGraph2D(ax2, 'TargAngl Hinge 3 [rad]', 'TargAngl Hinge 24 [rad]',mincoord = [closingAngl2[0], closingAngl1[0]], 
                maxcoord = [closingAngl2[-1], closingAngl1[-1]],  divisions = [divitheta2, divitheta1], buffer = [sep2, sep1])
else:
    NiceGraph2D(ax2, 'TargAngl Hinge 3 [rad]', 'TargAngl Hinge 24 [rad]',mincoord = [closingAngl1[0], closingAngl2[0]], 
                maxcoord = [closingAngl1[-1], closingAngl2[-1]],  divisions = [divitheta1, divitheta2], buffer = [sep1, sep2])

if inverted:
    cs2 = ax2.imshow(stableStateMat, extent=[theta2[0,0]-sep2,theta2[0,-1]+sep2,theta1[0,0]-sep1,theta1[-1,0]+sep1], 
                     cmap = cm.Set3, aspect = 'auto', origin = 'lower')
else:
    cs2 = ax2.imshow(stableStateMat, extent=[theta1[0,0]-sep1,theta1[-1,0]+sep1,theta2[0,0]-sep2,theta2[0,-1]+sep2], 
                     cmap = cm.Set3, aspect = 'auto', origin = 'lower')

ax2.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))
ax2.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))

cbar2 = plt.colorbar(cs2,  ax = ax2, fraction=0.05, pad=0.01)#, extend = 'max'
cbar2.set_ticks(np.linspace(0, max(inverse), max(inverse)+1))
cbar2.set_label('Stable State', fontsize = 15, color = '0.2')
cbar2.ax.tick_params(axis='y',colors='0.2')
cbar2.ax.tick_params(axis='x',colors='0.2')
cbar2.outline.set_edgecolor('0.2')

fig2.tight_layout()
fig2.show()
fig2.savefig(folder_name + '/StableStates.png', transparent = False)

#############################################################################
#Plot of real angles

realtheta1 = -dataAngles[2::IterPerSimul,2]/np.pi
realtheta1 = realtheta1[sortAngl[::-1]]
realtheta2 = -dataAngles[2::IterPerSimul,23]/np.pi
realtheta2 = realtheta2[sortAngl[::-1]]

fig3 = plt.figure(2,figsize=(cm2inch(24.1), cm2inch(20)))
ax3 = plt.subplot(111)
NiceGraph2D(ax3, 'Hinge 3 [rad]', 'Hinge 24 [rad]', mincoord = [min(realtheta1), min(realtheta2)], 
                maxcoord = [max(realtheta1), max(realtheta2)],  divisions = [divitheta2, divitheta1], buffer = [sep2, sep1])

cmap1, norm1 = from_levels_and_colors(np.linspace(0,maxEnergy,1000), cm.rainbow(np.linspace(0, 1, 1000-1)))
cmap1.set_over('r')
cs3 = ax3.scatter(realtheta1, realtheta2, c = totEnergysort, cmap = cmap1, vmax = maxEnergy, s = 360, marker = 's')

ax3.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))
ax3.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))

cbar3 = plt.colorbar(cs3, ax = ax3, fraction=0.05, pad=0.01, extend = 'max')
cbar3.set_ticks(np.linspace(0, maxEnergy, 5))
cbar3.set_label('Energy', fontsize = 15, color = '0.2')
cbar3.ax.tick_params(axis='y',colors='0.2')
cbar3.ax.tick_params(axis='x',colors='0.2')
cbar3.outline.set_edgecolor('0.2')

fig3.tight_layout()
fig3.show()
fig3.savefig(folder_name + '/RealAngles-Energy.png', transparent = False)