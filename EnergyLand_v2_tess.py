# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:04:51 2019

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

Folder_name = "Results/Tessellation4/25/2CFF/09-Jan-2020_2_1_"

allDesigns = pd.DataFrame()
# allKappasAnalysis = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    
tessellation = np.array(Folder_name.split('_')[1:-1]).astype(int)
numUC = np.prod(tessellation)
  
for subdir in os.listdir(Folder_name):
    if subdir == 'Images':
        continue    
        
    for subdir2 in os.listdir(Folder_name+'/'+subdir):
        folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
        
        ThisData, simLen = raa.ReadFile(folder_name)
        ThisData = raa.maskBadResults(ThisData, True)  
        
        ThisData.iloc[:,-numUC*4::] = raa.orderAnglesMat(ThisData.iloc[:,-numUC*4::].to_numpy(), 4, tessellation)
        ThisData['StableStates'] = raa.countStableStates(ThisData.iloc[:,-numUC*4::], 0.6, 'centroid')
        
        selData = raa.makeSelectionPerStSt(ThisData, simLen*0.0)
        
        allDesigns = allDesigns.append(selData)
        
    # kappasnumStSt = kappasStSt.groupby('kappa')['StableStates'].nunique()
    # kappasnumStSt = kappasnumStSt.reset_index()
    # kappasnumStSt['restang'] = np.ones(np.shape(kappasnumStSt)[0])*restang

    # allKappasAnalysis = allKappasAnalysis.append(kappasnumStSt)

allDesigns = allDesigns.reset_index(level=0, drop =True)

#%%
#### Get stable states from material
allDesigns['StableStateAll'] = raa.countStableStates(allDesigns.iloc[:,5:5+numUC*4].values, 1, 'centroid')
colormap = 'jet'

### Get stable states of individual vertices
posStSt = np.shape(allDesigns)[0]
allUCang = np.resize(allDesigns.iloc[:,8:8+numUC*4].values,(posStSt*numUC,4))
[allUCang,sym] = raa.orderAngles(allUCang, 4, posStSt*numUC)
matUCStSt = raa.countStableStates(allUCang, 0.7, 'centroid',True)
matUCStSt = raa.standarizeStableStates(matUCStSt, allUCang, onlysign = False)
matUCStSt = np.reshape(matUCStSt, (posStSt, numUC))
colormapUC = 'jet'

#### Get stable states from types of vertex
[allDesigns['StableStateFromUC'],allDesigns['StableStatefUCName']] = raa.extractStableStates(matUCStSt)
colormapfUC = 'Set3'

#%%
plt.close('all')   
    
##### Plotting the minFace of each stable state to make sure we are inside the constraint
plot.XYperZ(allDesigns, 0, r'$\kappa$', 7, r'$Area$', 1, -3, colormap, save = True, Folder_name = Folder_name, NameFig = 'AreaFaces')

#%%
plt.close('all')   

##### Plotting the Curvature of each stable state
ststcol = -2
plot.XYperZ(allDesigns, 0, r'$\kappa$', 6, r'$K$', 1, ststcol, colormapfUC, save = True, Folder_name = Folder_name, NameFig = 'Curvature')
plot.CreatLegend(allDesigns.iloc[:,ststcol+1], colormapfUC, save = True, Folder_name = Folder_name, NameFig = 'Curvature')

#%%
##### Plot first three angles of all vertices with different colors

##### Colors indicating the type of vertex
plot.Angles3D(allUCang, matUCStSt.flatten(), colormapUC)

# ##### Colors indicating the stable state of material
# plot.Angles3D(allUCang, np.tile(allDesigns['StableStateAll'].values,(numUC,1)).transpose().flatten(), colormap)

# ##### Colors indicating the stable state from types of vertex
# plot.Angles3D(allUCang, np.tile(allDesigns['StableStateFromUC'].values,(numUC,1)).transpose().flatten(), colormapfUC)

#%%
raa.SaveForPlot(allDesigns, Folder_name)

