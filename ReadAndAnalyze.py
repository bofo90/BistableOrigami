# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:48:23 2020

@author: iniguez
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matl
from matplotlib.colors import from_levels_and_colors
from mpl_toolkits import mplot3d
import configparser
import os.path
import copy

def ReadFile(folder_name):
    
    file_name1 = "/EnergyData.csv" 
    file_name2 = "/Hinges.csv"
    file_name3 = "/PosStad.csv"
    file_name4 = "/Angles.csv"
    file_metadata = "/metadata.txt"

    dataEnergy = pd.read_csv(folder_name+file_name1)
    simLen = np.size(dataEnergy['Hinge Number'])
    
    dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
    dataEnergy = dataEnergy.drop(['DiagonalEnergy','FaceEnergy','TargetAngleEnergy'], axis = 1)
    #Here needs to be an adjustment for old simulations
    dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']
    
    metadata = pd.DataFrame()
    [ang, kap] = ReadMetadata(folder_name+file_metadata)
    metadata['kappa'] = np.ones(simLen)*kap
    metadata['restang'] = np.ones(simLen)*ang
    
    dataCurv = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', dtype = np.float64)
    dataCurv = np.delete(dataCurv, 0, 1)
    dataEnergy['Curvature'] = np.mean(dataCurv,axis = 1)
    
    dataLen = np.loadtxt(folder_name+file_name3,skiprows=1, delimiter = ',', dtype = np.float64)
    dataLen = np.delete(dataLen, 0, 1)
    dataEnergy['minFace'] = np.min(dataLen,axis = 1)

    dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
    dataAngles = np.delete(dataAngles, 0, 1)

    allAngles = pd.DataFrame(dataAngles)
    dataEnergy = pd.concat([metadata, dataEnergy, allAngles], axis=1, sort=False)

    return dataEnergy, simLen

def maskBadResults(ThisData, printFlags = False):
    
    minFaceFlag = ThisData['minFace']<=0.11
    ThisData['Flags'][minFaceFlag] = -5
    
    if printFlags:
        print(ThisData['Flags'].value_counts())
    
    exfl = ThisData.Flags.values
    flagmask = (exfl !=1) & (exfl !=2) & (exfl !=0)
    # flagmask = ~flagmask.any(axis= 1)
    ThisData = ThisData.iloc[~flagmask,:]
    
    ThisData = ThisData.drop(['Flags'], axis = 1)
        
    return ThisData

def orderAngles(angles, ucsize, simulations):
    
    symetries = ["" for x in range(simulations)]
    finalAngles = np.zeros((simulations,np.size(angles,1)))
    for sim in np.arange(simulations):
        sym = ''
        ver = np.array(angles[sim,:])
        if np.sum(np.sign(np.around(ver, decimals = 0))) < 0: ### sort to have mayority positive angles
            # ver = ver[::-1]*(-1)
            ver = ver*(-1)
            sym = sym + 'm'
            
        rolnum = 4-np.argmin(ver)
        ver= np.roll(ver,rolnum)         ### sort from smallest to highest angles
        if rolnum != 4:
            sym = sym+'r'+ str(rolnum)
            
        finalAngles[sim,:] = ver
        sym = sym +'_'
        symetries[sim] = sym
        
    # return [angles, symetries]    
    return [finalAngles, symetries]

def orderAnglesMat(angles, ucsize, tessellation):
    
    nunitcells = np.prod(tessellation)
    simulations = np.size(angles,0)
    
    mirror = (np.sum(np.sign(np.around(angles, decimals = 4)),axis = 1) < 0)
    angles[mirror,:] = angles[mirror,:]*-1
    
    finalAngles = np.reshape(np.sum(np.reshape(angles,(simulations*nunitcells,4)),axis = 1),(simulations,nunitcells))
    
    orderedang = np.zeros(np.shape(angles))
    order = np.zeros([simulations, nunitcells]).astype(int)
    for i in np.arange(simulations):
        locorder = np.reshape(np.arange(nunitcells),tessellation)
        ststs = finalAngles[i,:]
        ststs = np.reshape(ststs, tessellation)
        rot = np.argmin([ststs[0,0], ststs[0,-1],ststs[-1,-1],ststs[-1,0]])
        ststs = np.rot90(ststs,rot)
        locorder = np.rot90(locorder,rot)
        mir = np.argsort([ststs[0,0], ststs[0,-1],ststs[-1,-1],ststs[-1,0]])
        if mir[1]> mir[2]:
            ststs = ststs.transpose()
            locorder = locorder.transpose()
        ststs = ststs.flatten()
        finalAngles[i,:] = ststs
        order[i,:] = locorder.flatten().astype(int)
        
        for j in np.arange(nunitcells):
            tempang = angles[i,order[i,j]*4:order[i,j]*4+4]
            tempang = np.roll(tempang, rot)
            if mir[1]> mir[2]:
                tempang = tempang[[1,0,3,2]]
            orderedang[i,j*4:j*4+4] = tempang
           
    return orderedang 
        
def countStableStates(finalAngles, distance, method, plot = False):
    
    import scipy.cluster.hierarchy as hierarch
    from scipy.spatial.distance import pdist
    
    if np.size(finalAngles,0)<=1:
        return [1]
                                       
    Z = hierarch.linkage(finalAngles, method)
    inverse = hierarch.fcluster(Z, distance, criterion='distance')
    c = hierarch.cophenet(Z, pdist(finalAngles))

    if plot:
        
        print('this is the cophenet of the hierarchical linkage', c[0])
        
        plt.figure(0,figsize=(10, 5))
        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('sample index')
        plt.ylabel('distance')
        hierarch.dendrogram(
            Z,
            truncate_mode='lastp',  # show only the last p merged clusters
            p=50,  # show only the last p merged clusters
            leaf_rotation=90.,  # rotates the x axis labels
            leaf_font_size=8.,  # font size for the x axis labels
            show_contracted=True,
            )
        plt.show()

    
    return inverse

def makeSelectionPerStSt(allData, simBound):
    
    
    kappasStSt = allData.groupby('StableStates').apply(lambda _df: _df.mean())
        
    kappasStSt['Hinge Number'] =  allData.groupby('StableStates')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 0)).to_numpy()
    
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('StableStates')[['Hinge Number']].count()
    more10percent = kappasStSt['amountStSt']>simBound
    kappasStSt = kappasStSt[more10percent]
   
    return kappasStSt


def getStableStatesMat(data, tessellation):
    
    numsim = np.size(data,0)
    stst = np.zeros(numsim)
    order = np.zeros([numsim, np.prod(tessellation)]).astype(int)
    
    for i in np.arange(numsim):
        locorder = np.reshape(np.arange(np.prod(tessellation)),tessellation)
        ststs = data[i,:]
        ststs = np.reshape(ststs, tessellation)
        rot = np.argmin([ststs[0,0], ststs[0,-1],ststs[-1,-1],ststs[-1,0]])
        ststs = np.rot90(ststs,rot)
        locorder = np.rot90(locorder,rot)
        mir = np.argsort([ststs[0,0], ststs[0,-1],ststs[-1,-1],ststs[-1,0]])
        if mir[1]> mir[2]:
            ststs = ststs.transpose()
            locorder = locorder.transpose()
        ststs = ststs.flatten().astype(int)
        colstst = np.array2string(ststs, separator='')
        stst[i] = colstst[1:-1]
        order[i,:] = locorder.flatten().astype(int)
        
    return stst, order

def standarizeStableStates(matUCStSt, allUCang, onlysign = False):
    
    standstst = np.round(np.array([[0,np.sqrt(2),0,np.sqrt(2)],[1,1,1,1],[-1,1,-1,1],[-1,1,1,1]])/2,1)
    
    ststUC = np.unique(matUCStSt)
    if onlysign:
        allUCang = np.sign(np.round(allUCang, 1))
    magangvec = np.sqrt(np.sum(allUCang*allUCang, 1))
    angvecunit = np.round(allUCang / magangvec[:,None],1)
    
    for i in np.arange(np.size(standstst,0)):
        ststloc = (angvecunit == standstst[i,:]).all(axis=1)
        oldStSt = np.unique(matUCStSt[ststloc])
        if (ststUC[oldStSt-1]<0).any():
            print("\n\nError: double standarizing!!!!! you need to check the standarization of stable states!!!\n\n")
        ststUC[oldStSt-1] = -i-1
    
    oldSS = np.unique(matUCStSt)
    [i, notdelSS, newSS] = np.unique(ststUC, return_index=True, return_inverse=True)
    newSS = newSS+1
    newStSt = np.zeros(np.shape(matUCStSt))
    for old,new in zip(oldSS, newSS):
        newStSt[matUCStSt == old] = new
        
    return newStSt.astype(int)

def extractStableStates(matUCStSt):
    
    simulations = np.size(matUCStSt,0)
    
    tempStSt = np.arange(simulations).astype(float)
    
    for sim in np.arange(simulations):
        unstst = np.unique(matUCStSt[sim,:])
        colstst = np.array2string(unstst, separator='')
        tempStSt[sim] = np.int(colstst[1:-1].replace(" ", ""))
        
    [i, newSS] = np.unique(tempStSt, return_inverse=True)
    newSS = newSS+1
    
    return newSS, tempStSt

def ReadMetadata(file):
    
    if os.path.isfile(file):
        metadata = configparser.RawConfigParser(allow_no_value=True)
        metadata.read(file)
    
        restang = float(metadata.get('options','restang'))
        kappa = float(metadata.get('options','Khinge'))
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

    return restang, kappa