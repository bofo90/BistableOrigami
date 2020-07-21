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

Folder_name = "Results/SingleVertex4/2CFF/24-Jun-2020_0.00_ 90.00_180.00_270.00_"

allDesignsSingleVer = pd.DataFrame()

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
        ThisData['StableStates'] = raa.countStableStatesDBSCAN(ThisData.iloc[:,-numvertex::], 0.1, 5)
        
        selData = raa.makeSelectionPerStSt(ThisData, simLen*0.0)
        clusData, clusFlags = raa.deleteNonClustered(selData,ThisFlags)
        
        allDesignsSingleVer = allDesignsSingleVer.append(clusData)
    
allDesignsSingleVer = allDesignsSingleVer.reset_index(level=0, drop =True)
allDesignsSingleVer = np.round(allDesignsSingleVer,8)

### Get stable states of vertices
allDesignsSingleVer['StableStateAll'] =raa.countStableStatesDistance(allDesignsSingleVer.iloc[:,8:8+numvertex].values, 1.5)


#%%
allDesigns = pd.DataFrame()
allCountMat = pd.DataFrame()
allMat = pd.DataFrame()
allEne = np.array([[0,0]])
allErrors = np.array([[0,0,0,0,0,0]])


for i in np.arange(4,5)[::-1]:

    Folder_name = "Results/Tessellation4/25/2CFF/19-Jun-2020_%d_%d_" %(i,i) #with no B.C.
    # Folder_name = "Results/Tessellation4/25/2CFF/07-Jul-2020_%d_%d_" %(i,i) #with no B.C.
    
    if not os.path.isdir(Folder_name):
        continue
    
    if not os.path.isdir(Folder_name + '/Images/'):
        os.mkdir(Folder_name + '/Images/')
           
    tessellation = np.array(Folder_name.split('_')[-3:-1]).astype(int)
    numUC = np.prod(tessellation)
    numvertex = np.int(Folder_name.split('/')[1][-1])
      
    for subdir in os.listdir(Folder_name):
        if subdir == 'Images':
            continue    
        # if subdir == 'RestAng_2.356': #'RestAng_1.571':
        #     continue
        # if subdir != 'RestAng_2.356': #'RestAng_0.785': #'RestAng_1.571': #
        #     continue
        if (subdir != 'RestAng_2.356') & (subdir != 'RestAng_0.785') & (subdir != 'RestAng_1.571'):
            continue
            
        for subdir2 in os.listdir(Folder_name+'/'+subdir):
            # if subdir2 !='kappa_0.03981':
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
            allAngles = np.resize(ThisAnglesMa,(simLen*numUC,numvertex))
            
            ### Different types of clustering the vertex' stable state
            # vertexStSt = raa.countStableStatesKmean(allAngles, 0.5, reduceFit = True)
            # vertexStSt = raa.countStableStatesDBSCAN(allAngles, 0.26, minpoints = 5, reduceFit = True)
            vertexStSt = raa.countStableStatesDistance(allAngles, 1.5)
            simStStMa = np.resize(vertexStSt, (simLen,numUC))
            
            maskPureMat, typePureMat, perPure, mat1Lines = raa.getPureMatComp(simStStMa, tessellation)
            ThisDataMa['StableStateMat'] = typePureMat
            ThisDataMa['Purity'] = perPure
            ThisDataMa['LineDefectMat1'] = np.argmax(mat1Lines, axis = 1)
            allErrors = np.concatenate((allErrors, mat1Lines), axis = 0)

            selDataMat = raa.makeSelectionPerStStMa(ThisDataMa)
            allMat = allMat.append(selDataMat)
            
            # plot.Angles3D(allAngles, vertexStSt, 'jet')
            simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure = raa.applyMask(maskPureMat, simStStMa, ThisDataMa, ThisEnergyMa, ThisAnglesMa, ThisCurvMa)
            
            selCountMat = raa.makeSelectionPureMat(ThisDataPure, ThisFlags, typePureMat)
            allCountMat = allCountMat.append(selCountMat)

            if ThisDataPure.empty:
                print("No pure material found")
                continue
            
            ### select just 500 simulations per parameter selection
            # selData = ThisDataMa.sample(500, random_state = 10)
            # allDesigns = allDesigns.append(selData)
            ### select all simulations per parameter selection
            allDesigns = allDesigns.append(ThisDataMa)
            
            
            ### save all energies
            # thisEne = np.concatenate(([np.ones(np.size(ThisEnergyMa.flatten()))*tessellation[0]], [ThisEnergyMa.flatten()]), axis = 0).T
            ### save pure energies
            # thisEne = np.concatenate(([np.ones(np.size(ThisEnergyPure.flatten()))*tessellation[0]], [ThisEnergyPure.flatten()]), axis = 0).T
            ### save non pure energies
            ThisEnergyNonPure = ThisEnergyMa[~maskPureMat,:]
            thisEne = np.concatenate(([np.ones(np.size(ThisEnergyNonPure.flatten()))*tessellation[0]], [ThisEnergyNonPure.flatten()]), axis = 0).T
            ### save center vertex energy
            # thisEneVert = ThisEnergyPure[:,np.int(np.floor(np.prod(tessellation)/2))]
            # thisEne = np.concatenate(([np.ones(np.size(thisEneVert))*tessellation[0]], [thisEneVert]), axis = 0).T
            allEne = np.concatenate((allEne, thisEne ), axis = 0)

        
allDesigns = allDesigns.reset_index(level=0, drop =True)
allCountMat = allCountMat.reset_index(level = 0, drop = True)
allMat = allMat.reset_index(level = 0, drop = True)
allEne = allEne[1:,:]

#%%
plt.close('all') 

##### plot appearance of material types with errors
plot.ColorbarPerZKappa(allCountMat,0, np.arange(9)+3, 1, save = True, Folder_name = Folder_name, NameFig = 'SimulationsConvergence')

#%%
plt.close('all')   
    
#### Plotting the Curvature, Energy and minFace of materials against kappa with color meaning max occurence of defect
thetas = np.unique(allDesigns.iloc[:,1])
for t in thetas:
    here = (allDesigns.iloc[:,1] == t).values
    hereSV = (allDesignsSingleVer.iloc[:,1] == np.round(t,8)).values
    plot.XKappawSingleVertPlot(allDesigns.iloc[here,:], allDesignsSingleVer.iloc[hereSV,:], 
                               0, 0, r'$\kappa$', 7, 6, r'$K_\mathregular{G}$', 11, 
                               save = True, Folder_name = Folder_name, NameFig = 'NeighvsCurvatureMat_ang%.2f_sel' %t)
    plot.XKappawSingleVertPlot(allDesigns.iloc[here,:], allDesignsSingleVer.iloc[hereSV,:],
                               0, 0, r'$\kappa$', 6, 5, r'$E_{norm}$', 11, 
                               save = True, Folder_name = Folder_name, NameFig = 'NeighvsEnergyMat_ang%.2f_sel' %t)
    plot.XKappawSingleVertPlot(allDesigns.iloc[here,:], allDesignsSingleVer.iloc[hereSV,:],
                               0, 0, r'$\kappa$', 8, 7, r'$min_Face$', 11, 
                               save = True, Folder_name = Folder_name, NameFig = 'NeighvsMinFaceMat_ang%.2f_sel' %t)

#%%
# allMat_copy = allMat.copy()
# allMat_copy['restang'] = allMat_copy['restang']*np.pi
# raa.SaveForPlotMat(allMat_copy, Folder_name[:42])

#### create concentric squares in matrix with increasing value
# x = 9
# a = np.zeros((x,x))
# for i in np.arange(x):
#     a[i,np.arange(i,x-i)] = i
#     a[x-1-i,np.arange(i,x-i)] = i
#     a[np.arange(i+1,x-1-i), i] = i
#     a[np.arange(i+1,x-1-i), x-1-i] = i