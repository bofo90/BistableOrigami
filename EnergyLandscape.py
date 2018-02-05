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

folder_name = "Results/truncated tetrahedron/active-set/energy/05-Feb-2018_Energylandscape_3to24\kh0.001_kta1.000_ke1.000_kf100.000"

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



#%%

divitheta1 = len(np.unique(closingAngl1))
divitheta2 = len(np.unique(closingAngl2))

theta1 = -closingAngl1.reshape((divitheta1,divitheta2))
theta2 = -closingAngl2.reshape((divitheta1,divitheta2))


fig1 = plt.figure(0,figsize=(cm2inch(35), cm2inch(20)))
ax1 = plt.subplot(111)
NiceGraph2D(ax1, 'Theta1 [rad]', 'Theta2 [rad]',mincoord = [-closingAngl1[0], -closingAngl2[0]], maxcoord = [-closingAngl1[-1], -closingAngl2[-1]],  divisions = [divitheta1, divitheta2])

cs1 = ax1.imshow(eTotalRel[1::2].reshape((divitheta1,divitheta2)).T, extent=[theta1[0,0],theta1[-1,0],theta2[0,0],theta2[0,-1]], cmap = cm.copper, aspect = 'auto')#theta1[:,0],theta2[0,:],
ax1.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))
ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g $\pi$'))

cbar = plt.colorbar(cs1,  ax = ax1, fraction=0.05, pad=0.01)#format="%d", extend = 'max'
cbar.set_ticks(np.linspace(0, max(eTotalRel[1::2]), 5))
cbar.set_label('Energy', fontsize = 15, color = '0.2')
cbar.ax.tick_params(axis='y',colors='0.2')
cbar.ax.tick_params(axis='x',colors='0.2')
cbar.outline.set_edgecolor('0.2')

fig1.tight_layout()
fig1.show()
fig1.savefig(folder_name + '/EnergyLand.png', transparent = False)








