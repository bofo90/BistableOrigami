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
allEne = np.array([[0,0]])


for i in np.arange(2,16)[::-1]:

    Folder_name = "Results/Tessellation4/25/2CFF/01-Apr-2020_%d_%d_" %(i,i) #with no B.C.
    # Folder_name = "Results/Tessellation4/25/2CFF/24-Apr-2020_%d_%d_" %(i,i) #with B.C.
    # Folder_name = "Results/Tessellation4/25/2CFF/03-Jun-2020_%d_%d_" %(i,i) #with new B.C.
    # Folder_name = "Results/Tessellation4/25/2CFF/08-May-2020_%d_%d_" %(i,i) #higher kappa
    # Folder_name = "Results/Tessellation4/25/2CFF/29-May-2020_%d_%d_" %(i,i) #higher kappa with P.B.C.
    
    if not os.path.isdir(Folder_name):
        continue
    
    if not os.path.isdir(Folder_name + '/Images/'):
        os.mkdir(Folder_name + '/Images/')
        
    # if i != 6:
    #     continue
        
    tessellation = np.array(Folder_name.split('_')[-3:-1]).astype(int)
    numUC = np.prod(tessellation)
    numvertex = np.int(Folder_name.split('/')[1][-1])
      
    for subdir in os.listdir(Folder_name):
        if subdir == 'Images':
            continue    
        if subdir == 'RestAng_1.571':
            continue
        # if subdir != 'RestAng_0.785': #'RestAng_2.356': #'RestAng_1.571': #
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
            # print('Converged simulations',  np.sum(ThisMask))
            
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
            maskPureMat, typePureMat, perPure, mat1Lines = raa.getPureMatConv(simStStMa, tessellation)
            ThisDataMa['StableStateMat'] = typePureMat
            ThisDataMa['Purity'] = perPure
            ThisDataMa['LineDefectMat1'] = mat1Lines

            selDataMat = raa.makeSelectionPerStStMa(ThisDataMa)
            allMat = allMat.append(selDataMat)
            
            # plot.Angles3D(allAngles, vertexStSt, 'jet')
            simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure = raa.applyMask(maskPureMat, simStStMa, ThisDataMa, ThisEnergyMa, ThisAnglesMa, ThisCurvMa)
            
            selCountMat = raa.makeSelectionPureMat(ThisDataPure, ThisFlags, typePureMat)
            allCountMat = allCountMat.append(selCountMat)

            if ThisDataPure.empty:
                print("No pure material found")
                continue
            
            # selData = ThisDataMa.sample(500, random_state = 10)
            # allDesigns = allDesigns.append(selData)
            allDesigns = allDesigns.append(ThisDataMa)
            
            
            # pos = [0,np.int(np.floor(tessellation[0]/2)),np.int(np.floor(np.prod(tessellation)/2))]
            # thisEne = np.concatenate(([np.ones(np.size(ThisEnergyMa.flatten()))*tessellation[0]], [ThisEnergyMa.flatten()]), axis = 0).T
            # thisEne = np.concatenate(([np.ones(np.size(ThisEnergyPure.flatten()))*tessellation[0]], [ThisEnergyPure.flatten()]), axis = 0).T
            ThisEnergyNonPure = ThisEnergyMa[~maskPureMat,:]
            thisEne = np.concatenate(([np.ones(np.size(ThisEnergyNonPure.flatten()))*tessellation[0]], [ThisEnergyNonPure.flatten()]), axis = 0).T
            # thisEneVert = ThisEnergyPure[:,np.int(np.floor(np.prod(tessellation)/2))]
            # thisEne = np.concatenate(([np.ones(np.size(thisEneVert))*tessellation[0]], [thisEneVert]), axis = 0).T
            allEne = np.concatenate((allEne, thisEne ), axis = 0)
            # selData = raa.makeSelectionVertMat(simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure, tessellation)
            
            # notOutLie = selData['StableStates'] != -1
            # allDesigns = allDesigns.append(selData.iloc[notOutLie.values,:])
        
allDesigns = allDesigns.reset_index(level=0, drop =True)
allCountMat = allCountMat.reset_index(level = 0, drop = True)
allMat = allMat.reset_index(level = 0, drop = True)
allEne = allEne[1:,:]
#%%
### Get stable states from material
# allDesAng = allDesigns.iloc[:,12:16].values
# ang_4D = raa.getAng4D(allDesAng)
# ang_4D_wpbc = np.concatenate((np.sin(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2),np.cos(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2)), axis = 1)
# allDesigns['StableStateAll'] = raa.countStableStatesDBSCAN(ang_4D_wpbc, 0.1,7)

# allDesAngOrd,sym = raa.orderAngles(allDesAng, numvertex, np.size(allDesAng,0))
# ang_4D = raa.getAng4D(allDesAngOrd)
# ang_4D_wpbc = np.concatenate((np.sin(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2),np.cos(ang_4D/[np.pi, np.pi, np.pi*2]*np.pi*2)), axis = 1)
# allDesigns['StableStateVert'] = raa.countStableStatesKmean(ang_4D_wpbc, 0.07)

# allDesigns['StableStateAll'] = allDesigns['StableStates'].astype(int)
# allMat['StableStateAll'] = allMat['StableStateMat'].astype(int)

colormap = 'jet'

# plot.Angles3D(allDesAng, allDesigns['StableStateAll'], colormap)
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
plt.close('all')   
    
#### Plotting the Curvature and Energy of materials against neighbours for restang
plot.XSizePlot(allDesigns, 2, r'$matSize$', 7, r'$K_\mathregular{G}$', 9, 11, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_sel')
plot.XSizePlot(allDesigns, 2, r'$matSize$', 6, r'$E_{norm}$', 9, 11, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_sel')


#%%
plt.close('all')   
    
##### Plotting the Curvature and Energy of materials against neighbours for restang
plot.violinPlot(allDesigns, 2, r'$matSize$', 7, r'$K_\mathregular{G}$', 9, 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_viol')
plot.violinPlot(allDesigns, 2, r'$matSize$', 6, r'$E_{norm}$', 9, 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_viol')

#%%
plt.close('all')  

plot.violin_scatter(allDesigns, 2, r'$matSize$', 7, r'$K_\mathregular{G}$', 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_comb')
plot.violin_scatter(allDesigns, 2, r'$matSize$', 6, r'$E_{norm}$', 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_comb')

#%%
plt.close('all')

# plot.XYperZ(allDesigns, 10, r'$Purity$', 6, r'$E_{norm}$', 2, 11, colormap, save = True, Folder_name = Folder_name, NameFig = 'PurityvsEnergyMat_size')
# plot.XYperZ(allDesigns, 10, r'$Purity$', 7, r'$K_\mathregular{G}$', 2, 11, colormap, save = True, Folder_name = Folder_name, NameFig = 'PurityvsCurvatureMat_size')
#%%
allMat_copy = allMat.copy()
allMat_copy['restang'] = allMat_copy['restang']*np.pi
raa.SaveForPlotMat(allMat_copy, Folder_name[:42])


x = 9
a = np.zeros((x,x))
for i in np.arange(x):
    a[i,np.arange(i,x-i)] = i
    a[x-1-i,np.arange(i,x-i)] = i
    a[np.arange(i+1,x-1-i), i] = i
    a[np.arange(i+1,x-1-i), x-1-i] = i