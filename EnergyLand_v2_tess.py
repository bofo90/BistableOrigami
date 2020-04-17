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

allDesigns = pd.DataFrame()
allFlags = pd.DataFrame()

for i in np.arange(2,5)[::-1]:

    Folder_name = "Results/Tessellation4/25/2CFF/01-Apr-2020_%d_%d_" %(i,i)
    
    if not os.path.isdir(Folder_name + '/Images/'):
        os.mkdir(Folder_name + '/Images/')
        
    tessellation = np.array(Folder_name.split('_')[1:-1]).astype(int)
    numUC = np.prod(tessellation)
    numvertex = np.int(Folder_name.split('/')[1][-1])
      
    for subdir in os.listdir(Folder_name):
        if subdir == 'Images':
            continue    
            
        for subdir2 in os.listdir(Folder_name+'/'+subdir):
            folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
            
            ThisEnergy, ThisAngles, ThisCurv, ThisData, simLen = raa.ReadFileMat(folder_name)
            ThisData, ThisMask = raa.maskBadResults(ThisData, returnMask = True)  
            
            ThisEnergy = ThisEnergy[ThisMask]
            ThisAngles = ThisAngles[ThisMask]
            ThisCurv = ThisCurv[ThisMask]
            
            simLen = np.size(ThisData,0)
            allAngles = np.resize(ThisAngles,(simLen*numUC,numvertex))
            allAngles,sym = raa.orderAngles(allAngles, numvertex, simLen*numUC)
            ThisAngles = np.resize(allAngles, (simLen,numUC*numvertex))
            vertexStSt = raa.countStableStatesDBSCAN(allAngles, 0.26)
            simStSt = np.resize(vertexStSt, (simLen,numUC))
            
            # plot.Angles3D(allAngles, vertexStSt, 'jet')
            
            
            selData = raa.makeSelectionVertMat(simStSt, ThisData, ThisEnergy, ThisAngles, ThisCurv, tessellation)
            notOutLie = selData['StableStates'] != -1
            
            allDesigns = allDesigns.append(selData.iloc[notOutLie.values,:])
        
allDesigns = allDesigns.reset_index(level=0, drop =True)

#%%
### Get stable states from material
allDesAng = allDesigns.iloc[:,12:16].values
ang_4D = raa.getAng4D(allDesAng)
ang_4D_wpbc = np.concatenate((np.sin(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2),np.cos(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2)), axis = 1)
allDesigns['StableStateAll'] = raa.countStableStatesDBSCAN(ang_4D_wpbc, 0.1,7)

colormap = 'Set2'

plot.Angles3D(allDesAng, allDesigns['StableStateAll'], colormap)

# ### Get stable states of individual vertices
# posStSt = np.shape(allDesigns)[0]
# allUCang = np.resize(allDesigns.iloc[:,8:8+numUC*4].values,(posStSt*numUC,4))
# [allUCang,sym] = raa.orderAngles(allUCang, 4, posStSt*numUC)
# matUCStSt = raa.countStableStates(allUCang, 0.7, 'centroid',True)
# matUCStSt = raa.standarizeStableStates(matUCStSt, allUCang, onlysign = False)
# matUCStSt = np.reshape(matUCStSt, (posStSt, numUC))
# colormapUC = 'jet'

# #### Get stable states from types of vertex
# [allDesigns['StableStateFromUC'],allDesigns['StableStatefUCName']] = raa.extractStableStates(matUCStSt)
# colormapfUC = 'Set3'

#%%
plt.close('all')   
    
##### Plotting the Curvature and Energy against tess for the different vertex
plot.XYperZ(allDesigns, 2, r'$x-rep$', 11, r'$K_\mathregular{G}$', -2, -1, colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')
plot.XYperZ(allDesigns, 2, r'$x-rep$', 10, r'$E_{norm}$', -2, -1, colormap, save = True, Folder_name = Folder_name, NameFig = 'Energy')
plot.CreateColorbar(allDesigns.iloc[:,-1], colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')

#%%
# raa.SaveForPlotMatt(allDesigns, Folder_name[:42])

