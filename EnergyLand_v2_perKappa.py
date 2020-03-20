# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:53:05 2020

@author: iniguez
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path
import ReadAndAnalyze as raa
import Plotting as plot


#%%
plt.close('all')

Folder_name = "Results/SingleVertex4/2CFF/20-Mar-2020_0.00_ 90.00_180.00_270.00_"

allDesigns = pd.DataFrame()
allFlags = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    
designang = np.array(Folder_name.split('_')[1:-1]).astype(float)
numvertex = np.int(Folder_name.split('/')[1][-1])
  
for subdir in os.listdir(Folder_name):
    if subdir == 'Images':
        continue      
    
    for subdir2 in os.listdir(Folder_name+'/'+subdir):
        folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
        
        ThisData, simLen = raa.ReadFile(folder_name)
        ThisData, ThisFlags = raa.maskBadResults(ThisData, returnStad = True) 
        
        simLen = np.size(ThisData,0)
        ThisData.iloc[:,-numvertex::],sym = raa.orderAngles(ThisData.iloc[:,-numvertex::].to_numpy(), numvertex, simLen)
        ThisData['StableStates'] = raa.countStableStates(ThisData.iloc[:,-numvertex::], 0.6, 'centroid')
        
        selData = raa.makeSelectionPerStSt(ThisData, simLen*0.0)
        
        allDesigns = allDesigns.append(selData)
        allFlags = allFlags.append(ThisFlags)
    
allDesigns = allDesigns.reset_index(level=0, drop =True)
allFlags = allFlags.reset_index(level=0, drop =True)
#%%

### Get stable states of vertices
allDesAng = allDesigns.iloc[:,8:8+numvertex].values
allDesigns['StableStateAll'] = raa.countStableStates(allDesAng, 0.7, 'centroid', True)
allDesigns['StableStateAll'] = raa.standarizeStableStates(allDesigns['StableStateAll'], allDesAng, onlysign = False)

colormap = 'Set2'
        
#%%
plt.close('all')    
    
##### Plotting the minFace of each stable state to make sure we are inside the constraint
plot.XYperZ(allDesigns, 1, r'$\theta_0$', 7, r'$Area$', 0, -1, colormap, save = True, Folder_name = Folder_name, NameFig = 'AreaFaces')


#%%
plt.close('all')   

##### Plotting the Curvature of each stable state
ststcol = -1
plot.XYperZ(allDesigns,  1, r'$\theta_0$', 6, r'$k_\mathregular{G}$', 0, ststcol, colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')
plot.CreateColorbar(allDesigns.iloc[:,ststcol], colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')

#%%
plt.close('all')   

##### Plotting the Curvature of each stable state
ststcol = -1
plot.XYperZ(allDesigns,1, r'$\theta_0$', 5, r'$E_{tot}$', 0, ststcol, colormap, save = True, Folder_name = Folder_name, NameFig = 'TotEnergyNorm')
plot.CreateColorbar(allDesigns.iloc[:,ststcol], colormap, save = True, Folder_name = Folder_name, NameFig = 'TotEnergyNorm')

plot.XmultYperZ(allDesigns, 1, r'$\theta_0$', [3,4], r'$E_{norm}$', 0, save = True, Folder_name = Folder_name, NameFig = 'Energies')
#%%
plt.close('all')  

##### Plot first three angles of all vertices with different colors representing the stable state
plot.Angles3D(allDesAng, allDesigns['StableStateAll'].values, colormap)

#%%
raa.SaveForPlot(allDesigns, Folder_name)

