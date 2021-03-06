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

# Folder_name = "Results/SingleVertex4/2CFF/06-Apr-2020_0.00_ 90.00_180.00_270.00_"
Folder_name = "Results/SingleVertex4/2CFF/19-Mar-2020_0.00_ 90.00_180.00_270.00_"
# Folder_name = "Results/SingleVertex4/2CFF/13-May-2020_d_0.00_ 90.00_180.00_270.00_"

allDesigns = pd.DataFrame()
allFlags = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    
designang = np.array(Folder_name.split('_')[-5:-1]).astype(float)
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
        # ThisData['StableStates'] = raa.countStableStatesKmean(ThisData.iloc[:,-numvertex::], 0.01)
        ThisData['StableStates'] = raa.countStableStatesDBSCAN(ThisData.iloc[:,-numvertex::], 0.1, 5)
        # ThisData['StableStates'] = raa.countStableStatesDistance(ThisData.iloc[:,-numvertex::], 0.01)
        
        selData = raa.makeSelectionPerStSt(ThisData, simLen*0.0)
        clusData, clusFlags = raa.deleteNonClustered(selData,ThisFlags)
        
        allDesigns = allDesigns.append(clusData)
        allFlags = allFlags.append(clusFlags)
    
allDesigns = allDesigns.reset_index(level=0, drop =True)
allFlags = allFlags.reset_index(level=0, drop =True)
#%%

### Get stable states of vertices
allDesAng = allDesigns.iloc[:,8:8+numvertex].values
# allDesigns['StableStateAll'] = raa.countstandardStableStates(allDesAng, onlysign = False)
# allDesigns['StableStateAll'] = raa.standarizeStableStates(allDesigns['StableStateAll'], allDesAng, onlysign = False)

# ang_4D = raa.getAng4D(allDesAng)
# ang_4D_wpbc = np.concatenate((np.sin(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2),np.cos(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2)), axis = 1)
# # allDesigns['StableStateAll'] = raa.countStableStatesKmean(ang_4D_wpbc, 0.07)

# allDesigns['StableStateAll'] = raa.countStableStatesDBSCAN(ang_4D_wpbc, 0.1,7)
# allDesigns['StableStateAll'] = raa.getFlatStates(allDesAng, allDesigns['StableStateAll'].values)

allDesigns['StableStateAll'] = raa.countStableStatesDistance(allDesAng, 1.5)

allFlagsRed , allDesignsRed = raa.makeSelectionClustered(allDesigns, allFlags, 50)

allDesAngRed = allDesignsRed.iloc[:,8:8+numvertex].values

allFlagsRed.iloc[:,1] /= 1000
allDesignsRed.iloc[:,13] /= 1000


colormap = 'Set2'

# plot.Angles3D(ang_4D_wpbc, allDesigns['StableStateAll'].values, colormap)
plot.Angles3D(allDesAng, allDesigns['StableStateAll'].values, colormap)
        
#%%
# plt.close('all')    
    
##### Plotting the minFace of each stable state to make sure we are inside the constraint
# plot.XYperZ(allDesigns, 1, r'$\Theta/\pi$', 7, r'$Area$', 0, -1, colormap, save = True, Folder_name = Folder_name, NameFig = 'AreaFaces')


#%%
plt.close('all')   

##### Plotting the Curvature of each stable state
ststcol = -1
# plot.CurvaturePaper(allDesigns,  1, r'$\Theta/\pi$', 6, r'$K_\mathregular{G}$', 0, ststcol, colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')
# plot.CreateColorbar(allDesigns.iloc[:,ststcol], colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')

#%%

# plt.close('all')   

# ##### Plotting the Curvature of each stable state
# ststcol = -1
# plot.XYperZ(allDesigns, 0, r'$\kappa$', 6, r'$K_\mathregular{G}$', 1, ststcol, colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')
# plot.CreateColorbar(allDesigns.iloc[:,ststcol], colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')

#%%

# plot.XYperZ(allDesigns,  1, r'$\theta_0/\pi$', 8, r'$\theta_1$', 0, -1, colormap, save = False, Folder_name = Folder_name, NameFig = 'Curvature')


#%%
# plt.close('all')   

# ##### Plotting the Curvature of each stable state
# ststcol = -1
# plot.TotEnergyperZ(allDesigns,1, r'$\theta_0/\pi$', 0, ststcol, colormap, save = True, Folder_name = Folder_name)
# plot.CreateColorbar(allDesigns.iloc[:,ststcol], colormap, save = True, Folder_name = Folder_name, NameFig = 'TotEnergyNorm')

# plot.XmultYperZ(allDesigns, 1, r'$\theta_0/\pi$', [3,4], r'$E_{norm}$', 0, save = True, Folder_name = Folder_name, NameFig = 'Energies')

#%%
plt.close('all')

###ploting scatter points of the presence of each stabel state
colormap = 'coolwarm'
plot.StableStatesCounturPlotCurvature(allDesignsRed,1, r'$\Theta/\pi$', 0, r'$\kappa$', -1, 6, r'$K_\mathregular{G}$', colormap, [-3,3], save = True, Folder_name = Folder_name, NameFig = 'StStAppearance')
colormap = 'Set2'

#%%
plt.close('all')

###ploting scatter points of the presence of each stabel state
colormap = 'viridis'
plot.StableStatesCounturPlotPaper(allDesignsRed,1, r'$\Theta/\pi$', 0, r'$\kappa$', -1, 5, r'$E$', colormap, [-3,0], save = True, Folder_name = Folder_name, NameFig = 'StStEnergy')
colormap = 'Set2'

#%%
plt.close('all')

###ploting scatter points of the presence of each stabel state
colormap = 'cividis_r'
plot.StableStatesCounturPlot(allFlagsRed,3, r'$\Theta/\pi$', 2, r'$\kappa$', 0, 1, r'$n/N_\mathregular{sim}$', colormap, [0,1], save = True, Folder_name = Folder_name, NameFig = 'FlagsAppearance')
plot.StableStatesCounturPlotPaper2(allDesignsRed,1, r'$\Theta/\pi$', 0, r'$\kappa$', -1, -2, r'$n/N_\mathregular{sim}$', colormap, [0,1], save = True, Folder_name = Folder_name, NameFig = 'StStNumSim')

colormap = 'Set2'

#%%
# plt.close('all')  

#### Plot first three angles of all vertices with different colors representing the stable state
plot.Angles3D(allDesAngRed, allDesignsRed['StableStateAll'].values, colormap)

#%%
allDes_copy = allDesignsRed.copy()
allDes_copy['restang'] = allDes_copy['restang']*np.pi
raa.SaveForPlot(allDes_copy, Folder_name)

