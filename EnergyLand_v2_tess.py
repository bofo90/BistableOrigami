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
        
    # if i != 15:
    #     continue
        
    tessellation = np.array(Folder_name.split('_')[-3:-1]).astype(int)
    numUC = np.prod(tessellation)
    numvertex = np.int(Folder_name.split('/')[1][-1])
      
    for subdir in os.listdir(Folder_name):
        if subdir == 'Images':
            continue    
        # if subdir == 'RestAng_1.571':
        #     continue
        # if subdir != 'RestAng_2.356': #'RestAng_1.571': #'RestAng_0.785': #
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
            vertexStSt = raa.countStableStatesDistance(allAngles, 10)
            simStStMa = np.resize(vertexStSt, (simLen,numUC))
            
            # maskPureMat, typePureMat = raa.getPureMat(simStStMa, tessellation)
            maskPureMat, typePureMat, perPure, mat1Lines, grainsize = raa.getPureMatComp(simStStMa, tessellation)
            ThisDataMa['StableStateMat'] = typePureMat
            ThisDataMa['Purity'] = perPure
            
            grainsizeCountMat = pd.DataFrame([grainsize.min(axis = 0)], columns = ['GSMat1', 'GSMat2', 'GSMat3']) 

            ThisDataDef = pd.concat([ThisDataMa.reset_index(level=0, drop =True),pd.DataFrame(mat1Lines), pd.DataFrame(grainsize, columns = ['GSMat1', 'GSMat2', 'GSMat3'])], axis=1, sort = True)
            selDataMat = raa.makeSelectionPerStStMa(ThisDataDef)
            allMat = allMat.append(selDataMat)
            
            # plot.Angles3D(allAngles, vertexStSt, 'jet')
            simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure = raa.applyMask(maskPureMat, simStStMa, ThisDataMa, ThisEnergyMa, ThisAnglesMa, ThisCurvMa)
            
            selCountMat = raa.makeSelectionPureMat(ThisDataPure, ThisFlags, typePureMat)
            allCountMat = allCountMat.append(pd.concat([selCountMat, grainsizeCountMat], axis=1, sort = True))

            # selData = ThisDataMa.sample(500, random_state = 10)
            # allDesigns = allDesigns.append(selData)
            allDesigns = allDesigns.append(ThisDataDef)

            if ThisDataPure.empty:
                print("No pure material found")
                continue
            
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
colormap = 'jet'

#%%
plt.close('all') 

plot.ColorbarPerZ(allCountMat,2, np.array([3,6,4,7,5,8,9,10]), 1, save = True, Folder_name = Folder_name, NameFig = 'SimulationsConvergence')

# #%%
# plt.close('all')   
    
# #### Plotting the Curvature and Energy of materials against neighbours for restang
# plot.XSizePlot(allDesigns, 2, r'$matSize$', 7, r'$K_\mathregular{G}$', 9, 10, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_sel')
# plot.XSizePlot(allDesigns, 2, r'$matSize$', 6, r'$E_{norm}$', 9, 10, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_sel')


# #%%
# plt.close('all')   
    
# ##### Plotting the Curvature and Energy of materials against neighbours for restang
# plot.violinPlot(allDesigns, 2, r'$matSize$', 7, r'$K_\mathregular{G}$', 9, 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_viol')
# plot.violinPlot(allDesigns, 2, r'$matSize$', 6, r'$E_{norm}$', 9, 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_viol')

# #%%
# plt.close('all')  

# plot.violin_scatter(allDesigns, 2, r'$matSize$', 7, r'$K_\mathregular{G}$', 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_comb')
# plot.violin_scatter(allDesigns, 2, r'$matSize$', 6, r'$E_{norm}$', 9, save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_comb')

# #%%
# plt.close('all')

# plot.XYperZ(allDesigns, 10, r'$Purity$', 6, r'$E_{norm}$', 2, 11, colormap, save = True, Folder_name = Folder_name, NameFig = 'PurityvsEnergyMat_size')
# plot.XYperZ(allDesigns, 10, r'$Purity$', 7, r'$K_\mathregular{G}$', 2, 11, colormap, save = True, Folder_name = Folder_name, NameFig = 'PurityvsCurvatureMat_size')

#%%
plt.close('all')

allMat = allMat.round(8)
thetas = np.unique(allMat.iloc[:,1])
for t in thetas:
    here = (allMat.iloc[:,1] == t).values
    plot.DefectsApperance(allMat.iloc[here,:], 2, r'$matSize$', save = True, Folder_name = Folder_name, NameFig = 'DefectsvsMatSize_ang%.2f_sel' %t)

#%%
plt.close('all')

#### Plotting kappa against num of simulations for all the different defects
allCountMat = allCountMat.round(8)
thetas = np.unique(allCountMat.iloc[:,1])
for t in thetas:
    here = (allCountMat.iloc[:,1] == t).values
    plot.GrainSize(allCountMat.iloc[here,:], 2, r'$matSize$', save = True, Folder_name = Folder_name, NameFig = 'GrainSizevsMatSize_ang%.2f' %t)

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