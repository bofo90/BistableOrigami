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

Folder_name = "Results/SingleVertex4/2CFF/02-Dec-2019_0.00_90.00_180.00_270.00_"

allDesigns = pd.DataFrame()
# allKappasAnalysis = pd.DataFrame()

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
        ThisData = raa.maskBadResults(ThisData) 
        
        simLen = np.size(ThisData,0)
        ThisData.iloc[:,-numvertex::],sym = raa.orderAngles(ThisData.iloc[:,-numvertex::].to_numpy(), numvertex, simLen)
        ThisData['StableStates'] = raa.countStableStates(ThisData.iloc[:,-numvertex::], 0.6, 'centroid')
        
        selData = raa.makeSelectionPerStSt(ThisData, simLen*0.0)
        
        allDesigns = allDesigns.append(selData)
        
    # kappasStSt['amountSym'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['Symmetry']].nunique())
    # kappasStSt['Symetries'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates').apply(lambda x: str(np.sort(x.Symmetry.unique()))))
    
    # kappasnumStSt = kappasStSt.groupby('kappa')['StableStates'].nunique()
    # kappasnumStSt = kappasnumStSt.reset_index()
    # kappasnumStSt['restang'] = np.ones(np.shape(kappasnumStSt)[0])*restang

    # allKappasAnalysis = allKappasAnalysis.append(kappasnumStSt)
    
allDesigns = allDesigns.reset_index(level=0, drop =True)

#%%

### Get stable states of vertices
allDesAng = allDesigns.iloc[:,8:8+numvertex].values
allDesigns['StableStateAll'] = raa.countStableStates(allDesAng, 0.7, 'centroid', True)
allDesigns['StableStateAll'] = raa.standarizeStableStates(allDesigns['StableStateAll'], allDesAng, onlysign = False)

colormap = 'Set2'

#### Mask results out from the specified range of kappas
# kappasrange = [10**-3,10**0]
# dropers = ~((allDesigns.kappa < kappasrange[0]) | (allDesigns.kappa > kappasrange[-1]))
# allDesigns = allDesigns[dropers]
# allDesAng = allDesAng[dropers]

#%%

# symstst = ["" for x in stst]

# for i in np.arange(np.size(stst,0)):
#     allsym = allDesigns[allDesigns['StableStateAll'] == stst[i]]['Symetries'].unique()
#     redsym = []
#     for sym in allsym:
#         redsym = np.append(redsym, sym[1:-1].split())
#     symstst[i] = ' '.join(np.unique(redsym))
#     print('Stable state',stst[i],'has symmetries:', symstst[i])
        
#%%
plt.close('all')    
    
##### Plotting the minFace of each stable state to make sure we are inside the constraint
plot.XYperZ(allDesigns, 0, r'$\kappa$', 7, r'$Area$', 1, -3, colormap, save = True, Folder_name = Folder_name, NameFig = 'AreaFaces')


#%%
plt.close('all')   

##### Plotting the Curvature of each stable state
ststcol = -1
plot.XYperZ(allDesigns, 0, r'$\kappa$', 6, r'$K_\mathregular{G}$', 1, ststcol, colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')
plot.CreateColorbar(allDesigns.iloc[:,ststcol], colormap, save = True, Folder_name = Folder_name, NameFig = 'Curvature')

#%%
plt.close('all')   

##### Plotting the Curvature of each stable state
ststcol = -1
plot.XYperZ(allDesigns, 0, r'$\kappa$', 5, r'$E_{tot}$', 1, ststcol, colormap, save = True, Folder_name = Folder_name, NameFig = 'TotEnergyNorm')
plot.CreateColorbar(allDesigns.iloc[:,ststcol], colormap, save = True, Folder_name = Folder_name, NameFig = 'TotEnergyNorm')

plot.XmultYperZ(allDesigns, 0, r'$\kappa$', [3,4], r'$E_{norm}$', 1, save = True, Folder_name = Folder_name, NameFig = 'Energies')
#%%
plt.close('all')  

##### Plot first three angles of all vertices with different colors representing the stable state
plot.Angles3D(allDesAng, allDesigns['StableStateAll'].values, colormap)

#%%
raa.SaveForPlot(allDesigns, Folder_name)

