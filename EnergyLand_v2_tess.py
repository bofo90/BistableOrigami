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
allCountMat = pd.DataFrame()
allMat = pd.DataFrame()

for i in np.arange(2,12)[::-1]:

    # Folder_name = "Results/Tessellation4/25/2CFF/01-Apr-2020_%d_%d_" %(i,i) #with no B.C.
    # Folder_name = "Results/Tessellation4/25/2CFF/24-Apr-2020_%d_%d_" %(i,i) #with B.C.
    # Folder_name = "Results/Tessellation4/25/2CFF/03-Jun-2020_%d_%d_" %(i,i) #with new B.C.
    Folder_name = "Results/Tessellation4/25/2CFF/08-May-2020_%d_%d_" %(i,i) #higher kappa
    # Folder_name = "Results/Tessellation4/25/2CFF/29-May-2020_%d_%d_" %(i,i) #higher kappa with P.B.C.
    
    if not os.path.isdir(Folder_name + '/Images/'):
        os.mkdir(Folder_name + '/Images/')
        
    tessellation = np.array(Folder_name.split('_')[-3:-1]).astype(int)
    numUC = np.prod(tessellation)
    numvertex = np.int(Folder_name.split('/')[1][-1])
      
    for subdir in os.listdir(Folder_name):
        if subdir == 'Images':
            continue    
        # if subdir != 'RestAng_1.571': #'RestAng_2.356': #'RestAng_0.785': #
        #     continue
            
        for subdir2 in os.listdir(Folder_name+'/'+subdir):
            # if subdir2 =='kappa_0.31623':
            #     continue
            folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
            print('Analysing: ' + folder_name)
            
            ThisEnergyOR, ThisAnglesOR, ThisCurvOR, ThisDataOR, simLen = raa.ReadFileMat(folder_name)
            ThisDataOR["TotalEnergy"] = np.mean(ThisEnergyOR, axis = 1)
            ThisDataMa, ThisMask, ThisFlags = raa.maskBadResults(ThisDataOR, returnMask = True, returnStad = True)  
            ThisFlags['tes'] = tessellation[0]
            print('Converged simulations',  np.sum(ThisMask))
            
            ThisEnergyMa = ThisEnergyOR[ThisMask]
            ThisAnglesMa = ThisAnglesOR[ThisMask]
            ThisCurvMa = ThisCurvOR[ThisMask]
            
            simLen = np.size(ThisDataMa,0)
            # allAngles = raa.selectSpecialVertex(ThisAngles, simLen, tessellation, numvertex)
            # allAngles,sym = raa.orderAngles(allAngles, numvertex, simLen*3)
            # ThisAngles = np.resize(allAngles, (simLen,3*numvertex))
            # vertexStSt = raa.countStableStatesDBSCAN(allAngles, 0.26)
            # simStSt = np.resize(vertexStSt, (simLen,3))
            
            allAngles = np.resize(ThisAnglesMa,(simLen*numUC,numvertex))
            # allAngles,sym = raa.orderAngles(allAngles, numvertex, simLen*numUC)
            # ThisAngles = np.resize(allAngles, (simLen,numUC*numvertex))
            
            # vertexStSt = raa.countStableStatesKmean(allAngles, 0.5, reduceFit = True)
            # vertexStSt = raa.countStableStatesDBSCAN(allAngles, 0.26, minpoints = 5, reduceFit = True)
            vertexStSt = raa.countStableStatesDistance(allAngles, 1.5)
            simStStMa = np.resize(vertexStSt, (simLen,numUC))
            
            # maskPureMat, typePureMat = raa.getPureMat(simStStMa, tessellation)
            maskPureMat, typePureMat = raa.getPureMatConv(simStStMa, tessellation)
            ThisDataMa['StableStateMat'] = typePureMat

            selDataMat = raa.makeSelectionPerStStMa(ThisDataMa)
            allMat = allMat.append(selDataMat)
            
            # plot.Angles3D(allAngles, vertexStSt, 'jet')
            simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure = raa.applyMask(maskPureMat, simStStMa, ThisDataMa, ThisEnergyMa, ThisAnglesMa, ThisCurvMa)
            
            selCountMat = raa.makeSelectionPureMat(ThisDataPure, ThisFlags, typePureMat)
            allCountMat = allCountMat.append(selCountMat)

            if ThisDataPure.empty:
                print("No pure material found")
                continue
            
            selData = raa.makeSelectionVertMat(simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure, tessellation)
            
            notOutLie = selData['StableStates'] != -1
            allDesigns = allDesigns.append(selData.iloc[notOutLie.values,:])
        
allDesigns = allDesigns.reset_index(level=0, drop =True)
allCountMat = allCountMat.reset_index(level = 0, drop = True)
allMat = allMat.reset_index(level = 0, drop = True)

#%%
### Get stable states from material
allDesAng = allDesigns.iloc[:,12:16].values
# ang_4D = raa.getAng4D(allDesAng)
# ang_4D_wpbc = np.concatenate((np.sin(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2),np.cos(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2)), axis = 1)
# allDesigns['StableStateAll'] = raa.countStableStatesDBSCAN(ang_4D_wpbc, 0.1,7)

# allDesAngOrd,sym = raa.orderAngles(allDesAng, numvertex, np.size(allDesAng,0))
# ang_4D = raa.getAng4D(allDesAngOrd)
# ang_4D_wpbc = np.concatenate((np.sin(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2),np.cos(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2)), axis = 1)
# allDesigns['StableStateVert'] = raa.countStableStatesKmean(ang_4D_wpbc, 0.07)

allDesigns['StableStateAll'] = allDesigns['StableStates'].astype(int)
allMat['StableStateAll'] = allMat['StableStateMat'].astype(int)

colormap = 'jet'

plot.Angles3D(allDesAng, allDesigns['StableStateAll'], colormap)
# plot.Angles3D(allDesAngOrd, allDesigns['StableStateVert'], colormap)

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

plot.ColorbarPerZ(allCountMat,2, np.arange(9)+3, 1, save = True, Folder_name = Folder_name, NameFig = 'SimulationsConvergence')

#%%
# plt.close('all')   
    
##### Plotting the Curvature and Energy of materials against neighbours for restang
plot.XYperZwError(allMat, 2, r'$matSize$', 7, r'$K_\mathregular{G}$', 1, -1, colormap, 11, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_ang')
plot.XYperZwError(allMat, 2, r'$matSize$', 6, r'$E_{norm}$', 1, -1, colormap, 10, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_ang')
plot.CreateColorbar(allMat.iloc[:,-1], colormap, save = True, Folder_name = Folder_name, NameFig = 'StableStates')

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
# plt.close('all')  

# colormap = 'jet'
# plot.XYperZline(allDesigns, 2, r'$matSize$', -4, r'$Amount pureStSt$', -3, 9, colormap, save = True, Folder_name = Folder_name, NameFig = 'NeighvsSim_sv')
# plot.XYperZ(allDesigns, 2, r'$matSize$', 11, r'$E_{norm}$', -3, 9, colormap, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergy_sv')
# plot.CreateColorbar(allDesigns.iloc[:,9], colormap, save = True, Folder_name = Folder_name, NameFig = 'StableStatesAll')
# colormap = 'Set2'
#%%
# plt.close('all')   
    
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
allMat_copy = allMat.copy()
allMat_copy['restang'] = allMat_copy['restang']*np.pi
raa.SaveForPlotMat(allMat_copy, Folder_name[:42])

