# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 12:13:51 2020

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


for i in np.arange(1,16)[::-1]:

    Folder_name = "Results/Tessellation4/25/2CFF/06-Jul-2020_1_%d_" %i 

    
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
        # if subdir == 'RestAng_1.571':
        #     continue
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
            ThisFlags['tes'] = tessellation[1]
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
#             maskPureMat, typePureMat, perPure, mat1Lines = raa.getPureMatLine(simStStMa, tessellation)
#             ThisDataMa['StableStateMat'] = typePureMat
#             ThisDataMa['Purity'] = perPure
#             ThisDataMa['LineDefectMat1'] = mat1Lines

#             selDataMat = raa.makeSelectionPerStStMa(ThisDataMa)
#             allMat = allMat.append(selDataMat)
            
#             # plot.Angles3D(allAngles, vertexStSt, 'jet')
#             simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure = raa.applyMask(maskPureMat, simStStMa, ThisDataMa, ThisEnergyMa, ThisAnglesMa, ThisCurvMa)
            
#             selCountMat = raa.makeSelectionPureMat(ThisDataPure, ThisFlags, typePureMat)
#             allCountMat = allCountMat.append(selCountMat)

#             if ThisDataPure.empty:
#                 print("No pure material found")
#                 continue
            
#             # selData = ThisDataMa.sample(500, random_state = 10)
#             # allDesigns = allDesigns.append(selData)
#             allDesigns = allDesigns.append(ThisDataMa)
            
            
#             # pos = [0,np.int(np.floor(tessellation[0]/2)),np.int(np.floor(np.prod(tessellation)/2))]
#             # thisEne = np.concatenate(([np.ones(np.size(ThisEnergyMa.flatten()))*tessellation[0]], [ThisEnergyMa.flatten()]), axis = 0).T
#             # thisEne = np.concatenate(([np.ones(np.size(ThisEnergyPure.flatten()))*tessellation[0]], [ThisEnergyPure.flatten()]), axis = 0).T
#             ThisEnergyNonPure = ThisEnergyMa[~maskPureMat,:]
#             thisEne = np.concatenate(([np.ones(np.size(ThisEnergyNonPure.flatten()))*tessellation[0]], [ThisEnergyNonPure.flatten()]), axis = 0).T
#             # thisEneVert = ThisEnergyPure[:,np.int(np.floor(np.prod(tessellation)/2))]
#             # thisEne = np.concatenate(([np.ones(np.size(thisEneVert))*tessellation[0]], [thisEneVert]), axis = 0).T
#             allEne = np.concatenate((allEne, thisEne ), axis = 0)
#             # selData = raa.makeSelectionVertMat(simStStPure, ThisDataPure, ThisEnergyPure, ThisAnglesPure, ThisCurvPure, tessellation)
            
#             # notOutLie = selData['StableStates'] != -1
#             # allDesigns = allDesigns.append(selData.iloc[notOutLie.values,:])
        
# allDesigns = allDesigns.reset_index(level=0, drop =True)
# allCountMat = allCountMat.reset_index(level = 0, drop = True)
# allMat = allMat.reset_index(level = 0, drop = True)
# allEne = allEne[1:,:]
