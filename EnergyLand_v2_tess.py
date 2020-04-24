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

for i in np.arange(1,15)[::-1]:

    Folder_name = "Results/Tessellation4/25/2CFF/01-Apr-2020_%d_%d_" %(i,i)
    
    if not os.path.isdir(Folder_name + '/Images/'):
        os.mkdir(Folder_name + '/Images/')
        
    tessellation = np.array(Folder_name.split('_')[1:-1]).astype(int)
    numUC = np.prod(tessellation)
    numvertex = np.int(Folder_name.split('/')[1][-1])
      
    for subdir in os.listdir(Folder_name):
        if subdir == 'Images':
            continue    
        # if subdir != 'RestAng_0.785': #'RestAng_1.571': #'RestAng_2.356': #
        #     continue
            
        for subdir2 in os.listdir(Folder_name+'/'+subdir):
            if subdir2 =='kappa_0.31623':
                continue
            folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
            print('Analysing: ' + folder_name)
            
            ThisEnergy, ThisAngles, ThisCurv, ThisData, simLen = raa.ReadFileMat(folder_name)
            ThisData["TotalEnergy"] = np.mean(ThisEnergy, axis = 1)
            ThisData, ThisMask = raa.maskBadResults(ThisData, returnMask = True)  
            
            ThisEnergy = ThisEnergy[ThisMask]
            ThisAngles = ThisAngles[ThisMask]
            ThisCurv = ThisCurv[ThisMask]
            
            simLen = np.size(ThisData,0)
            # allAngles = raa.selectSpecialVertex(ThisAngles, simLen, tessellation, numvertex)
            # allAngles,sym = raa.orderAngles(allAngles, numvertex, simLen*3)
            # ThisAngles = np.resize(allAngles, (simLen,3*numvertex))
            # vertexStSt = raa.countStableStatesDBSCAN(allAngles, 0.26)
            # simStSt = np.resize(vertexStSt, (simLen,3))
            
            allAngles = np.resize(ThisAngles,(simLen*numUC,numvertex))
            # allAngles,sym = raa.orderAngles(allAngles, numvertex, simLen*numUC)
            ThisAngles2 = np.resize(allAngles, (simLen,numUC*numvertex))
            vertexStSt = raa.countStableStatesDBSCAN(allAngles, 0.26, reduceFit = True)
            simStSt2 = np.resize(vertexStSt, (simLen,numUC))
            
            simStSt, ThisData, ThisEnergy, ThisAngles1, ThisCurv = raa.getPureMat(simStSt2, ThisData, ThisEnergy, ThisAngles2, ThisCurv, tessellation)
            
            # plot.Angles3D(allAngles, vertexStSt, 'jet')
            
            if ThisData.empty:
                print("No pure material found")
                continue
            selData = raa.makeSelectionVertMat(simStSt, ThisData, ThisEnergy, ThisAngles1, ThisCurv, tessellation)
            
            notOutLie = selData['StableStates'] != -1
            allDesigns = allDesigns.append(selData.iloc[notOutLie.values,:])
            
            # allDesigns = allDesigns.append(selData)
        
allDesigns = allDesigns.reset_index(level=0, drop =True)

#%%
### Get stable states from material
allDesAng = allDesigns.iloc[:,12:16].values
# allDesAng,sym = raa.orderAngles(allDesAng, numvertex, np.size(allDesAng,0))
ang_4D = raa.getAng4D(allDesAng)
ang_4D_wpbc = np.concatenate((np.sin(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2),np.cos(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2)), axis = 1)
allDesigns['StableStateAll'] = raa.countStableStatesDBSCAN(ang_4D_wpbc, 0.1,7)
# allDesigns['StableStateAll'] = raa.countStableStatesKmean(ang_4D_wpbc, 0.07)

colormap = 'Set2'

plot.Angles3D(ang_4D_wpbc, allDesigns['StableStateAll'], colormap)
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
# plt.close('all')   
    
# ##### Plotting the Curvature and Energy against neighbours for restang
# plot.XYperZwError(allDesigns, 2, r'$matSize$', 11, r'$K_\mathregular{G}$', 1, -1, colormap, 17, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvature_ang')
# plot.XYperZwError(allDesigns, 2, r'$matSize$', 10, r'$E_{norm}$', 1, -1, colormap, 16, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergy_ang')
# plot.CreateColorbar(allDesigns.iloc[:,-1], colormap, save = True, Folder_name = Folder_name, NameFig = 'StableStates')

#%%
# colormap = 'jet'
# plot.XYperZwDoubleError(allDesigns, 6, r'$E_{norm}^\mathregular{mat}$', 10, r'$E_{norm}$', 1, 2, colormap, 18, 16, save = True, Folder_name = Folder_name, NameFig = 'EnergyMatvsVert_ang')
# plot.CreateColorbar(allDesigns.iloc[:,2], colormap, save = True, Folder_name = Folder_name, NameFig = 'MaterialSize')
# colormap = 'Set2'
#%%
plt.close('all')  

colormap = 'jet'
plot.XYperZline(allDesigns, 2, r'$matSize$', -4, r'$Amount pureStSt$', -3, -1, colormap, save = True, Folder_name = Folder_name, NameFig = 'NeighvsSim_sv')
plot.XYperZline(allDesigns, 2, r'$matSize$', 10, r'$E_{norm}$', -3, -1, colormap, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergy_sv')
plot.CreateColorbar(allDesigns.iloc[:,-1], colormap, save = True, Folder_name = Folder_name, NameFig = 'StableStatesAll')
colormap = 'Set2'
#%%
plt.close('all')   
    
##### Plotting the Curvature and Energy against tess for the different vertex
# plot.XYperZwError(allDesigns, -2, r'$neigh$', 11, r'$K_\mathregular{G}$', -3, -1, colormap, 17, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvature_sv')
# plot.XYperZwError(allDesigns, -2, r'$neigh$', 10, r'$E_{norm}$', -3, -1, colormap, 16, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergy_sv')

#%%
# plt.close('all')   
    
# colormap = 'jet'
# ##### Plotting the Curvature and Energy against tess for the different vertex
# plot.XYperZwDoubleError(allDesigns, 7, r'$K_\mathregular{G}^\mathregular{mat}$', 11, r'$K_\mathregular{G}$', -3, -2, colormap, 19, 17,save = True, Folder_name = Folder_name, NameFig = 'CurvatureMatvsVert_sv')
# plot.XYperZwDoubleError(allDesigns, 6, r'$E_{norm}^\mathregular{mat}$', 10, r'$E_{norm}$', -3, -2, colormap, 18, 16, save = True, Folder_name = Folder_name, NameFig = 'EnergyMatvsVert_sv')
# plot.CreateColorbar(allDesigns.iloc[:,-2], colormap, save = True, Folder_name = Folder_name, NameFig = 'Neigbourgs')
# colormap = 'Set2'

#%%

# plot.XYperZwError(allDesigns, 2, r'$matSize$', 11, r'$K_\mathregular{G}$', -1, -3, colormap, 17, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvature_stst')
# plot.XYperZwError(allDesigns, 2, r'$matSize$', 10, r'$E_{norm}$', -1, -3, colormap, 16, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergy_stst')
# plot.CreateColorbar(allDesigns.iloc[:,-3], colormap, save = True, Folder_name = Folder_name, NameFig = 'VertexType')


#%%
# plt.close('all')   
    
# colormap = 'jet'
# ##### Plotting the Curvature and Energy against tess for the different vertex
# plot.XYperZwDoubleError(allDesigns, 7, r'$K_\mathregular{G}^\mathregular{mat}$', 11, r'$K_\mathregular{G}$', -3, -2, colormap, 19, 17,save = True, Folder_name = Folder_name, NameFig = 'CurvatureMatvsVert_sv')
# plot.XYperZwDoubleError(allDesigns, 6, r'$E_{norm}^\mathregular{mat}$', 10, r'$E_{norm}$', -3, -2, colormap, 18, 16, save = True, Folder_name = Folder_name, NameFig = 'EnergyMatvsVert_sv')
# colormap = 'Set2'


#%%
# raa.SaveForPlotMatt(allDesigns, Folder_name[:42])

