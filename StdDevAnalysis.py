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

possibleStDv = np.arange(0.05,0.55,0.05)

allDesigns = pd.DataFrame()
allFlags = pd.DataFrame()

for stddev in possibleStDv:

    Folder_name = "Results/SingleVertex4/2CFF/28-Feb-2020_%.2f_StDev" %stddev
    
    if not os.path.isdir(Folder_name + '/Images/'):
        os.mkdir(Folder_name + '/Images/')
        
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
            
            selData['StdDev'] = np.ones((np.size(selData,0),1))*stddev
            ThisFlags['StdDev'] = np.ones((np.size(ThisFlags,0),1))*stddev
            
            allDesigns = allDesigns.append(selData)
            allFlags = allFlags.append(ThisFlags)
    
allDesigns = allDesigns.reset_index(level=0, drop =True)
allFlags = allFlags.reset_index(level=0, drop =True)
#%%
cutkappa = 0.1
cutang = 1.571#2.356#0.785
redDesignsMask = (np.round(allDesigns['kappa'],3)==cutkappa) & (np.round(allDesigns['restang'],3)==cutang)
redFlagsMask = (np.round(allFlags['kappa'],3)==cutkappa) & (np.round(allFlags['restang'],3)==cutang)
redDesign = allDesigns.iloc[redDesignsMask.to_numpy(),:]
redFlags = allFlags.iloc[redFlagsMask.to_numpy(),:]

#%%
redDesAng = redDesign.iloc[:,8:8+numvertex].values
redDesign['StableStateAll'] = raa.countStableStates(redDesAng, 0.7, 'centroid')
redDesign['StableStateAll'] = raa.standarizeStableStates(redDesign['StableStateAll'], redDesAng, onlysign = False)
redDesign.iloc[(redDesign['StableStateAll'] > 6).to_numpy(),-1] = 7

colormap = 'Set2'
#%%
##### Plot first three angles of all vertices with different colors representing the stable state
plot.Angles3D(redDesAng, redDesign['StableStateAll'].values, colormap)

#%%
##### Plotting the standard deviation against the stable states to see when does it appears
plot.XYperZ(redDesign, -2, r'$StdDev$', -1, r'$StableState$', 1, -1, colormap, save = False, Folder_name = Folder_name, NameFig = 'AreaFaces')


#%%
###### Plotting the Std Dev against the amount of stable state per stable state and against the amount of flags per flag
# plot.XYperZ(redDesign, -2, r'$StdDev$', -3, r'$AmountStSt$', -1, -1, colormap, save = False, Folder_name = Folder_name, NameFig = 'AreaFaces')
# plot.XYperZ(redFlags, -1, r'$StdDev$', 1, r'$AmountFLags$', 0, 0, 'jet', save = False, Folder_name = Folder_name, NameFig = 'AreaFaces')

# plot.ConvSim(redDesign,redFlags, 'jet')
plot.ConvSim(redDesign,redFlags, 'Set3')






