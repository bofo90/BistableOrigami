# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:48:23 2020

@author: iniguez
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import configparser
import os.path
from datetime import datetime

def ReadFileMat(folder_name):
    
    file_name1 = "/EnergyData.csv" 
    file_name2 = "/Hinges.csv"
    file_name3 = "/PosStad.csv"
    file_name4 = "/Angles.csv"
    file_name5 = "/Stretch.csv"
    file_metadata = "/metadata.txt"
    
    numhin, numedg, restang, numUC = ReadMetadataMat(folder_name+file_metadata)
    
    dataParam = pd.read_csv(folder_name+file_name1)
    simLen = np.size(dataParam['Hinge Number'])   
    
    metadata = pd.DataFrame()
    [ang, kap] = getRAandK(folder_name)
    metadata['kappa'] = np.ones(simLen)*kap
    if oldSample(folder_name):
        metadata['kappa'] = metadata['kappa']/4
    metadata['restang'] = np.ones(simLen)*ang/np.pi
    metadata['tes'] = np.int(np.sqrt(numUC))
        
    dataParam['TotalEnergy'] = dataParam['EdgeEnergy']+dataParam['DiagonalEnergy']+dataParam['HingeEnergy']
    dataParam = dataParam.drop(['DiagonalEnergy','FaceEnergy','TargetAngleEnergy'], axis = 1)
    dataParam['TotalEnergy'] = dataParam['TotalEnergy'] /(kap*numhin+numedg)
    dataParam['HingeEnergy'] = dataParam['HingeEnergy']/kap/numhin
    dataParam['EdgeEnergy'] = dataParam['EdgeEnergy']/numedg
    
    dataCurv = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', dtype = np.float64)
    dataCurv = np.delete(dataCurv, 0, 1)
    dataParam['Curvature'] = np.mean(dataCurv,axis = 1)
    
    dataArea = np.loadtxt(folder_name+file_name3,skiprows=1, delimiter = ',', dtype = np.float64)
    dataArea = np.delete(dataArea, 0, 1)
    dataParam['minFace'] = np.min(dataArea,axis = 1)

    dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
    dataAngles = np.delete(dataAngles, 0, 1)
    
    dataLen = np.loadtxt(folder_name+file_name5,skiprows=1, delimiter = ',', dtype = np.float64)
    dataLen = np.delete(dataLen,0,1)
    
    dataLen = np.reshape(dataLen, (simLen,numUC,8))
    dataAng = np.reshape(dataAngles, (simLen,numUC,4))
    dataEnergy = calculateEnergy(dataLen, dataAng, kap, restang)

    # allAngles = pd.DataFrame(dataAngles)
    # dataEnergy = pd.concat([metadata, dataEnergy, allAngles], axis=1, sort=False)
    dataParam = pd.concat([metadata, dataParam], axis=1, sort=False)

    return dataEnergy, dataAngles, dataCurv, dataParam, simLen

def calculateEnergy(Len, Ang, kap, restang):
    
    LenEnergy = 0.5*Len**2
    AngEnergy = kap/(restang**4)*(Ang**2-restang**2)**2
    
    Energy = np.sum(LenEnergy, axis = 2)+np.sum(AngEnergy, axis = 2)
    
    return Energy

def ReadFile(folder_name):
    
    file_name1 = "/EnergyData.csv" 
    file_name2 = "/Hinges.csv"
    file_name3 = "/PosStad.csv"
    file_name4 = "/Angles.csv"
    file_metadata = "/metadata.txt"
    
    dataEnergy = pd.read_csv(folder_name+file_name1)
    simLen = np.size(dataEnergy['Hinge Number'])   
    
    metadata = pd.DataFrame()
    [ang, kap] = getRAandK(folder_name)
    metadata['kappa'] = np.ones(simLen)*kap
    if oldSample(folder_name):
        metadata['kappa'] = metadata['kappa']/4
    metadata['restang'] = np.ones(simLen)*ang/np.pi   
    
    numhin, numedg = ReadMetadata(folder_name+file_metadata)
        
    dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
    dataEnergy = dataEnergy.drop(['DiagonalEnergy','FaceEnergy','TargetAngleEnergy'], axis = 1)
    dataEnergy['TotalEnergy'] = dataEnergy['TotalEnergy'] /(kap*numhin+numedg)
    dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']/kap/numhin
    dataEnergy['EdgeEnergy'] = dataEnergy['EdgeEnergy']/numedg
    
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

def oldSample(folder_name):
    
    splitfolder = folder_name.split('/')
    
    if splitfolder[1][:-1] == 'SingleVertex':
        date = splitfolder[3].split('_')[0]
    if splitfolder[1][:-1] == 'Tessellation':
        date = splitfolder[4].split('_')[0]
    
    simdate = datetime.strptime(date, '%d-%b-%Y').date()
    changedate = datetime(2020,1,8).date()
    
    if simdate <= changedate:
        # print('\nMaking change of kappa due to old results\n\n')
        return True
    else:
        return False
        

def maskBadResults(ThisData, printFlags = False, returnStad = False, returnMask = False):
    
    minFaceFlag = (ThisData['minFace']<=0.105).values
    
    if printFlags:
        print(ThisData['Flags'].value_counts())
         
    exfl = ThisData.Flags.values
    flagmask = (exfl !=1) & (exfl !=2) & (exfl !=5) #& (exfl !=0)
    # flagmask = ~flagmask.any(axis= 1)
    
    onlyminArea = minFaceFlag & ~flagmask  
    ThisData.loc[onlyminArea,'Flags'] = -5
    
    MaskedData = ThisData.iloc[~(flagmask | minFaceFlag),:]
    MaskedData = MaskedData.drop(['Flags'], axis = 1)
    # MaskedData = ThisData
    
    ThisData.iloc[~(flagmask | minFaceFlag), 5] = 1
    flagmask2 = (flagmask | minFaceFlag) & ~onlyminArea
    # flagmask2 = (flagmask & ~minFaceFlag)
    ThisData.iloc[flagmask2, 5] = -1
    
    flagStad = ThisData['Flags'].value_counts()
    flagStad = flagStad.reset_index(level=0)
    flagStad = flagStad.rename(columns={"index": "Flags", "Flags": "amountFlags"})
    flagStad['kappa'] = np.ones((np.size(flagStad,0),1))*ThisData.iloc[0,0]
    flagStad['restang'] = np.ones((np.size(flagStad,0),1))*ThisData.iloc[0,1]
    
    if returnMask:
        return MaskedData, ~(flagmask | minFaceFlag)
    
    if returnStad:
        return MaskedData, flagStad
    else:
        return MaskedData

def orderAngles(angles, ucsize, simulations):
    
    symetries = ["" for x in range(simulations)]
    finalAngles = np.zeros((simulations,np.size(angles,1)))
    for sim in np.arange(simulations):
        sym = ''
        ver = np.array(angles[sim,:])
        if np.sum(np.around(ver, decimals = 1)) < 0: ### sort to have mayority positive angles
            # ver = ver[::-1]*(-1)
            ver = ver*(-1)
            sym = sym + 'm'
            
        rolnum = np.size(angles,1)-np.argmin(ver)
        ver= np.roll(ver,rolnum)         ### sort from smallest to highest angles
        if rolnum != np.size(angles,1):
            sym = sym+'r'+ str(rolnum)
            
        mirr2 = np.argsort(ver)
        if mirr2[1]>mirr2[3]:
            ver[[1,3]]= ver[[3,1]]
            sym = sym+'am'
            
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

def countStableStatesKmean(finalAngles, dis):
    
    from sklearn.cluster import KMeans
    
    N = np.size(finalAngles,0)
    
    if N <= 1:
        return [1]
    
    for i in np.arange(500)+1:        
        kmeans = KMeans(n_clusters=i, max_iter=1000, verbose = False, random_state = 0 )
        kmeans.fit_transform(finalAngles)
        if np.sqrt(kmeans.inertia_/(N-1)) < dis: ### the stad. dev. of angles in one cluster should not be more than dis
            break

    if i == 500:
        print('The clustering has reached its maximum\n')
    inverse = kmeans.labels_
    return inverse

def countStableStatesDBSCAN(finalAngles, dis = 0.05, minpoints = 20):
    
    from sklearn.cluster import DBSCAN
    
    N = np.size(finalAngles,0)
    
    if N <= 1:
        return [1]
    
    db = DBSCAN(eps=dis, min_samples=minpoints).fit(finalAngles)
    
    inverse = db.labels_
    # if (inverse == -1).any():
    #     inverse = inverse+1
    
    return inverse

def makeSelectionPerStSt(allData, simBound):
    
    
    kappasStSt = allData.groupby('StableStates').apply(lambda _df: _df.mean())
        
    kappasStSt['Hinge Number'] =  allData.groupby('StableStates')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 0)).to_numpy()
    
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('StableStates')[['Hinge Number']].count()
    more10percent = kappasStSt['amountStSt']>simBound
    kappasStSt = kappasStSt[more10percent]
   
    return kappasStSt

def makeSelectionVertMat(matstst, allData, Energy, Angles, Curv, tess):
    
    allSpec = pd.DataFrame()
    
    pos = [0,np.int(np.floor(tess[0]/2)),np.int(np.floor(tess[0]/2)),np.int(np.floor(np.prod(tess)/2))]
    
    for j in np.arange(3):
        i = pos[j]
        specData = SelectPerVertex(matstst[:,i], allData, Energy[:,i], Angles[:,i*4:(i+1)*4], Curv[:,i])
        specData["VertexType"] = j
        
        allSpec = allSpec.append(specData)
   
    return allSpec

def SelectPerVertex(VertStSt, allData_or, VertEnergy, VertAngles, VertCurv):
    
    allData = allData_or.copy()
    allData = allData.reset_index(level=0, drop =True)
    
    allData['StableStates'] = VertStSt
    allData['VertEnergy'] = VertEnergy
    allData['VertCurv'] = VertCurv
    
    VertAngles = pd.DataFrame(VertAngles)
    allData = pd.concat([allData, VertAngles], axis=1, sort=False)
    
    
    selStSt = allData.groupby('StableStates').apply(lambda _df: _df.mean())
    
    error = allData.groupby('StableStates').apply(lambda _df: _df.std())
    selStSt['StdVertEnergy'] = error['VertEnergy'].values
    selStSt['StdVertCurv'] = error['VertCurv'].values
        
    selStSt['Hinge Number'] =  allData.groupby('StableStates')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 0)).to_numpy()
    
    ####Only select stable states that have more than 10% appearance in each simulation
    selStSt['amountStSt'] = allData.groupby('StableStates')[['Hinge Number']].count()
    # more10percent = kappasStSt['amountStSt']>simBound
    # kappasStSt = kappasStSt[more10percent]
    
    
    return selStSt

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
    
    standstst = np.round(np.array([[0,np.sqrt(2),0,np.sqrt(2)],[1,1,1,1],[-1,1,-1,1],[-1,1,1,1],[-1,-1,1,1]])/2,1)
    
    ststUC = np.unique(matUCStSt)
    if onlysign:
        allUCang = np.sign(np.round(allUCang, 1))
    magangvec = np.sqrt(np.sum(allUCang*allUCang, 1))
    angvecunit = np.round(allUCang / magangvec[:,None],1)
    
    for i in np.arange(np.size(standstst,0)):
        projection = np.sqrt(np.absolute(np.sum(angvecunit*standstst[i,:],1)))
        ststloc = (projection > 0.97) & (magangvec > 0.4)
        oldStSt = np.unique(matUCStSt[ststloc])
        if (ststUC[oldStSt-1]<0).any():
            print("\nError: double standarizing!!!!! you need to check the standarization of the typical stable states!!!\n")
        ststUC[oldStSt-1] = -i-1
        
    ##### get stable states around the the flat configuration
    # ststloc = magangvec <= 0.4
    # oldStSt = np.unique(matUCStSt[ststloc])
    # if (ststUC[oldStSt-1]<0).any():
    #     print("\nError: double standarizing!!!!! you need to check the standarization of the flat stable states!!!\n")
    # ststUC[oldStSt-1] = -np.size(standstst,0)-1
        
    
    oldSS = np.unique(matUCStSt)
    [i, notdelSS, newSS] = np.unique(ststUC, return_index=True, return_inverse=True)
    newSS = newSS+1
    newStSt = np.zeros(np.shape(matUCStSt))
    for old,new in zip(oldSS, newSS):
        newStSt[matUCStSt == old] = new
        
    return newStSt.astype(int)

def countstandardStableStates(allUCang, onlysign = False):
    
    standstst = np.round(np.array([[0,np.sqrt(2),0,np.sqrt(2)],[1,1,1,1],[-1,1,-1,1],[-1,1,1,1],[-1,-1,1,1]])/2,1)
    
    newStSt = np.zeros(np.size(allUCang,0))
    
    if onlysign:
        allUCang = np.sign(np.round(allUCang, 1))
    magangvec = np.sqrt(np.sum(allUCang*allUCang, 1))
    angvecunit = np.round(allUCang / magangvec[:,None],3)
    
    # for i in np.arange(np.size(standstst,0)):
    #     projection = np.sqrt(np.absolute(np.sum(angvecunit*standstst[i,:],1)))
    #     ststloc = (projection > 0.9612)
    #     if (newStSt[ststloc] != 0).any():
    #         print("\nError: double standarizing!!!!! you need to check the standarization of the typical stable states!!!\n")
    #     newStSt[ststloc] = i+1
    
    for i in np.arange(np.size(allUCang,0)):
        projection = np.sqrt(np.absolute(np.sum(angvecunit[i,:]*standstst,1)))
        ststloc = np.argmax(projection)
        if (projection[ststloc] > 0.95):
            newStSt[i] = ststloc+1
        
    ##### get stable states around the the flat configuration
    # ststloc = magangvec <= 0.4
    # if (newStSt[ststloc] != 0).any():
    #     print("\nError: double standarizing!!!!! you need to check the standarization of the flat stable states!!!\n")
    # newStSt[ststloc] = i+2
        
    [i, notdelSS, newSS] = np.unique(newStSt, return_index=True, return_inverse=True)
    newSSnames = np.arange(np.size(i))+1
        
    return newSSnames[newSS].astype(int)  

def getAng4D(angles):
    
    angles = np.roll(angles, 0,axis = 1)
    
    magangvec = np.sqrt(np.sum(angles*angles, 1))
    angvecunit = angles / magangvec[:,None]
    
    ang_4D = np.zeros((np.size(angles,0),3))
    
    for i in np.arange(2):
        ang_4D[:,i] = np.arctan(np.sqrt(np.sum(angvecunit[:,i+1:]**2,1))/angvecunit[:,i])+np.pi/2
        
    ang_4D[:,-1] = 2*np.arctan(angvecunit[:,-1]/np.sqrt(np.sum(angvecunit[:,-2:]**2,1))+angvecunit[:,-2])+np.pi
    
    return ang_4D

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
    
        hing = float(metadata.get('extUnitCell','hinges'))
        edge = float(metadata.get('extUnitCell','edges'))
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

    return hing, edge

def ReadMetadataMat(file):
    
    if os.path.isfile(file):
        metadata = configparser.RawConfigParser(allow_no_value=True)
        metadata.read(file)
    
        hing = float(metadata.get('extUnitCell','hinges'))
        edge = float(metadata.get('extUnitCell','edges'))
        restang = float(metadata.get('options','restang'))
        x = float(metadata.get('options','xrep'))
        y = float(metadata.get('options','yrep'))
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

    return hing, edge, restang, np.int(x*y)

def getRAandK(folder_name):
    
    splited = folder_name.split('/')
    
    kappa_name = splited[-2]
    restang_name = splited[-3]
    
    kap = np.float(kappa_name.split('_')[-1])
    ang = np.float(restang_name.split('_')[-1])
    
    return ang, kap

def SaveForPlot(s, Folder_name):
    
    des = s.copy()
    
    if oldSample(Folder_name):
        des['kappa'] = des['kappa']*4
    
    des[['kappa','Hinge Number','StableStateAll','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforAllImages.csv', index = False)
    des.groupby('StableStateAll').apply(lambda df: df.sample(1, random_state = 0))[['kappa','Hinge Number','StableStateAll','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforStableStatesImages.csv', index = False)
       
    return

def SaveForPlotMatt(allDesigns, folder):
    
    SaveForPlot(thisDesign, thisfolder)
    
    return
