# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:35:56 2020

@author: iniguez
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path
import ReadAndAnalyze as raa
import Plotting as plot

plt.close('all')

possibleStDv = np.linspace(0.05,0.5,10)

allDesigns = pd.DataFrame()

for stddev in possibleStDv:

    Folder_name = "Results/SingleVertex4/2CFF/24-Feb-2020_%.2f_StDev" %stddev
    
    if not os.path.isdir(Folder_name + '/Images/'):
        os.mkdir(Folder_name + '/Images/')
        
    numvertex = np.int(Folder_name.split('/')[1][-1])
      
    for subdir in os.listdir(Folder_name):
        if subdir == 'Images':
            continue      
        
        for subdir2 in os.listdir(Folder_name+'/'+subdir):
            folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
            
            ThisData, simLen = raa.ReadFile(folder_name)
            ThisData = raa.maskBadResults(ThisData, True) 
            
            simLen = np.size(ThisData,0)
            ThisData.iloc[:,-numvertex::],sym = raa.orderAngles(ThisData.iloc[:,-numvertex::].to_numpy(), numvertex, simLen)
            ThisData['StableStates'] = raa.countStableStates(ThisData.iloc[:,-numvertex::], 0.6, 'centroid')
            
            selData = raa.makeSelectionPerStSt(ThisData, simLen*0.0)
            
            selData['StdDev'] = np.ones((np.size(selData,0),1))*stddev
            
            allDesigns = allDesigns.append(selData)
    
allDesigns = allDesigns.reset_index(level=0, drop =True)
#%%
allDesAng = allDesigns.iloc[:,8:8+numvertex].values
allDesigns['StableStateAll'] = raa.countStableStates(allDesAng, 0.7, 'centroid', True)
allDesigns['StableStateAll'] = raa.standarizeStableStates(allDesigns['StableStateAll'], allDesAng, onlysign = False)

colormap = 'Set2'
#%%
##### Plot first three angles of all vertices with different colors representing the stable state
plot.Angles3D(allDesAng, allDesigns['StableStateAll'].values, colormap)

#%%
##### Plotting the minFace of each stable state to make sure we are inside the constraint
plot.XYperZ(allDesigns, -2, r'$StdDev$', -3, r'$AmountStSt$', -1, -1, colormap, save = False, Folder_name = Folder_name, NameFig = 'AreaFaces')









