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
    Energy = Energy/(kap*4+8)
    
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
    
    #### asign -5 to converged results with min area
    onlyminArea = minFaceFlag & ~flagmask  
    ThisData.loc[onlyminArea,'Flags'] = -5
    
    #### mask the results that either don't have a convered flag or have min area
    MaskedData = ThisData.iloc[~(flagmask | minFaceFlag),:]
    MaskedData = MaskedData.drop(['Flags'], axis = 1)
    # MaskedData = ThisData
    
    #### asign flag 1 to convered simulations and -1 to non-converged ones
    ThisData.loc[~(flagmask | minFaceFlag), 'Flags'] = 1
    flagmask2 = (flagmask | minFaceFlag) & ~onlyminArea
    # flagmask2 = (flagmask & ~minFaceFlag)
    ThisData.loc[flagmask2, 'Flags'] = -1
    
    flagStad = ThisData['Flags'].value_counts()
    flagStad = flagStad.reset_index(level=0)
    flagStad = flagStad.rename(columns={"index": "Flags", "Flags": "amountFlags"})
    flagStad['kappa'] = np.ones((np.size(flagStad,0),1))*ThisData.iloc[0,0]
    flagStad['restang'] = np.ones((np.size(flagStad,0),1))*ThisData.iloc[0,1]
    
    if returnMask & ~returnStad:
        return MaskedData, ~(flagmask | minFaceFlag)
    if returnStad & ~returnMask:
        return MaskedData, flagStad
    if returnStad & returnMask:
        return MaskedData, ~(flagmask | minFaceFlag), flagStad
    else:
        return MaskedData
    
def selectSpecialVertex(ThisAngles, simLen, tess, numvertex):
    
    angles = np.reshape(ThisAngles,(simLen,np.prod(tess), numvertex))
    pos = [0,np.int(np.floor(tess[0]/2)),np.int(np.floor(np.prod(tess)/2))]
    
    angVertex = np.reshape(angles[:,pos,:],(3*simLen,numvertex))
    
    return angVertex

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

def countStableStatesKmean(finalAngles, dis, reduceFit = False):
    
    from sklearn.cluster import KMeans
    
    N = np.size(finalAngles,0)
    
    if N <= 1:
        return [1]
    
    numFit = 1000
    if reduceFit and (N>numFit):
        np.random.seed(0)
        fittingsamples = np.random.randint(np.size(finalAngles,0), size = numFit)
        fitsamples = finalAngles[fittingsamples,:]
    else:
        fitsamples = finalAngles
    
    for i in np.arange(500)+1: 
        print('clusters: ',i)
        kmeans = KMeans(n_clusters=i, max_iter=1000, verbose = False, random_state = 0 )
        kmeans.fit(fitsamples)
        if np.sqrt(kmeans.inertia_/(N-1)) < dis: ### the stad. dev. of angles in one cluster should not be more than dis
            break
        
    if i == 500:
        print('The clustering has reached its maximum\n')
        
    
    inverse = kmeans.predict(finalAngles)
    return inverse

def countStableStatesDBSCAN(finalAngles, dis = 0.05, minpoints = 20, reduceFit = False):
    
    from sklearn.cluster import DBSCAN
    
    N = np.size(finalAngles,0)
    
    if N <= 1:
        return [1]
    
    numFit = 5000
    
    if reduceFit and (N>numFit):
        np.random.seed(0)
        fittingsamples = np.random.randint(np.size(finalAngles,0), size = numFit)
        
        db = DBSCAN(eps=dis, min_samples=minpoints).fit(finalAngles[fittingsamples,:])
        
        inverse = dbscan_predict(db, finalAngles)
        
    else:
        db = DBSCAN(eps=dis, min_samples=minpoints).fit(finalAngles)
        
        inverse = db.labels_
    
    return inverse

def countStableStatesDistance(finalAngles, dis = 0.05):
    
    startVect = np.zeros((14,4))
    startVect[10,:] = [np.pi, 0, np.pi, 0]
    startVect[11,:] = [-np.pi, 0, -np.pi, 0]
    startVect[12,:] = [0, np.pi, 0, np.pi]
    startVect[13,:] = [0, -np.pi, 0, -np.pi]
    
    direction = np.array([[1,1,1,1], [1,-1,1,-1],[1,1,1,-1],[1,-1,1,1],[1,1,-1,1],[-1,1,1,1],[1,1,-1,-1],[1,-1,-1,1],[1,0,1,0],[0,1,0,1],[0,-1,0,1],[0,-1,0,1],[-1,0,1,0],[-1,0,1,0]]).astype(float)
    direction = (direction.T/np.linalg.norm(direction,axis = 1)).T
    
    ba = direction
    magn_ba = np.linalg.norm(ba, axis = 1)
    
    N = np.size(finalAngles,0)
    # if N <= 1:
    #     return [1]
    
    distances = np.zeros((N,14))
    for i in np.arange(14):
        pa = finalAngles - startVect[i,:]
        t = np.sum(pa*ba[i,:],axis = 1)/magn_ba[i]**2
        distances[:,i] = np.linalg.norm(pa-np.array([t]).T*ba[i,:], axis = 1)
        
    inverse = np.argmin(distances,axis = 1)
    
    no_state = np.min(distances,axis = 1) > dis
    inverse[no_state] = -1
    
    magn = np.sqrt(np.sum(finalAngles**2,axis = 1))
    inverse[magn<0.1] = 14
    
    return inverse

def dbscan_predict(model, X):

    nr_samples = X.shape[0]

    y_new = np.ones(shape=nr_samples, dtype=int) * -1

    for i in range(nr_samples):
        diff = model.components_ - X[i, :]  # NumPy broadcasting

        dist = np.linalg.norm(diff, axis=1)  # Euclidean distance

        shortest_dist_idx = np.argmin(dist)

        if dist[shortest_dist_idx] < model.eps:
            y_new[i] = model.labels_[model.core_sample_indices_[shortest_dist_idx]]

    return y_new

def getFlatStates(angles, inverse):
    
    magn = np.sqrt(np.sum(angles**2,axis = 1))
    
    maxStSt = np.max(inverse)
    
    inverse[magn<0.1] = maxStSt+1
    
    return inverse

def getPureMat(simStSt, tessellation):
    
    purity = np.zeros(np.size(simStSt,0))
    candidate = np.zeros(np.size(simStSt,0))
    matName = np.zeros(np.size(simStSt,0))
    
    stst = np.unique(simStSt)
    if stst[0] == -1:
        stst = np.delete(stst,0)
        
    #checkerboard pattern (dome-saddle)
    x = np.zeros(tessellation, dtype = int) 
    x[1::2, ::2] = 1
    x[::2, 1::2] = 1
    # x[::2, 1::2] = 2
    # x[1::2, 1::2] = 3
    x = x.flatten()
    
    #lines (miura)
    yh = np.zeros(tessellation, dtype = int) 
    yv = np.zeros(tessellation, dtype = int) 
    yh[::2,:] = 1
    yv[:,::2] = 1
    yh = yh.flatten()
    yv = yv.flatten()
    
    #all same (fold)
    z = np.zeros(tessellation, dtype = int) 
    z = z.flatten()
    
    for i in stst:
        here = (simStSt[:,x == 0] == i).all(axis = 1)
        orhere = (simStSt[:,x == 1] == i).all(axis = 1)
        defhere = np.logical_xor(here, orhere)
        # orhere2 = (simStSt[:,x == 2] == i).all(axis = 1)
        # orhere3 = (simStSt[:,x == 3] == i).all(axis = 1)
        # defhere = np.logical_xor(np.logical_xor(here, orhere),np.logical_xor(orhere2, orhere3))
        candidate[defhere] = candidate[defhere] + 2*0.5
       
        here = (simStSt[:,yh == 0] == i).all(axis = 1)
        orhere = (simStSt[:,yh == 1] == i).all(axis = 1)
        defhere = np.logical_xor(here, orhere)
        candidate[defhere] = candidate[defhere] + 3*0.5
        
        here = (simStSt[:,yv == 0] == i).all(axis = 1)
        orhere = (simStSt[:,yv == 1] == i).all(axis = 1)
        defhere = np.logical_xor(here, orhere)
        candidate[defhere] = candidate[defhere] + 5*0.5
        
        here = (simStSt[:,z == 0] == i).all(axis = 1)
        candidate[here] = candidate[here] + 7*1
    
    match = candidate != 0
    matchName = matName[match]
    matchName[(candidate[match] %2) == 0] = 1
    matchName[(candidate[match] %3) == 0] = 2
    matchName[(candidate[match] %5) == 0] = 3
    matchName[(candidate[match] %7) == 0] = 4
    matName[match] = matchName 
    
    purity = matName != 0
    
    print('Not pure materials ', np.sum(~purity))
    
    return purity, matName

def getPureMatConv(simStSt, tessellation):
    
    from scipy import signal
    
    numUCmat = np.prod(tessellation-1)
    simulations = np.size(simStSt,0)
    simStSt = simStSt.reshape((simulations,*tessellation))
    matType = np.zeros((simulations,4))
    
    #convolution matriz to identify checkerboard pattern
    convMat1 = np.array([[1,2],[-2,-1]])
    convMat2 = np.array([[1,2],[-1,-2]])
    convMat3 = np.array([[1,-1],[2,-2]])
    convMat4 = np.array([[1,1],[1,1]])
    
    for i in np.arange(simulations):
        sumConMat = signal.convolve2d(simStSt[i,:,:], convMat4, mode = 'valid')
        
        #check for checkerboar pattern with dome saddle
        convMat = signal.convolve2d(simStSt[i,:,:], convMat1, mode = 'valid')
        matType[i,0] = np.sum((convMat == 0) & (sumConMat == 2))/numUCmat
        
        #check for vert lines pattern with miura from elastic origami
        convMat = signal.convolve2d(simStSt[i,:,:], convMat2, mode = 'valid')
        matType[i,1] = np.sum((convMat == 0) & (sumConMat == 18))/numUCmat
        #check for vert lines pattern with miura from rigid origami
        matType[i,3] = np.sum((convMat == 0) & (sumConMat == 42))/numUCmat
        
        #check for horz lines pattern with miura from elastic origami
        convMat = signal.convolve2d(simStSt[i,:,:], convMat3, mode = 'valid')
        matType[i,1] += np.sum((convMat == 0) & (sumConMat == 10))/numUCmat
        #check for horz lines pattern with miura from rigid origami
        matType[i,3] += np.sum((convMat == 0) & (sumConMat == 50))/numUCmat
        
        #check for all same with fold from rigid origami
        matType[i,2] = np.sum(sumConMat == 32)/numUCmat
        matType[i,2] += np.sum(sumConMat == 36)/numUCmat
   
    purity = np.max(matType, axis = 1) == 1
    matName = np.argmax(matType, axis = 1)+1
    matName[~purity] += 4
    
    nomat = np.max(matType, axis = 1) == 0
    matName[nomat] = 0
    
    print('Not pure materials ', np.sum(~purity))
    
    return purity, matName
    
def applyMask(purity, simStSt, ThisData, ThisEnergy, ThisAngles, ThisCurv):
    
    simStSt = simStSt[purity,:]
    ThisData = ThisData.iloc[purity,:]
    ThisEnergy = ThisEnergy[purity,:]
    ThisAngles = ThisAngles[purity,:]
    ThisCurv = ThisCurv[purity,:]
    
    return simStSt, ThisData, ThisEnergy, ThisAngles, ThisCurv

def makeSelectionPerStSt(allData, simBound):
    
    
    kappasStSt = allData.groupby('StableStates').apply(lambda _df: _df.mean())
        
    kappasStSt['Hinge Number'] =  allData.groupby('StableStates')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 0)).to_numpy()
    
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('StableStates')[['Hinge Number']].count()
    more10percent = kappasStSt['amountStSt']>simBound
    kappasStSt = kappasStSt[more10percent]
   
    return kappasStSt

def makeSelectionPerStStMa(allData, simBound = 0):
    
    
    kappasStSt = allData.groupby('StableStateMat').apply(lambda _df: _df.mean())
    
    error = allData.groupby('StableStateMat').apply(lambda _df: _df.std())
    kappasStSt['StdMatEnergy'] = error['TotalEnergy'].values
    kappasStSt['StdMatCurv'] = error['Curvature'].values
        
    kappasStSt['Hinge Number'] =  allData.groupby('StableStateMat')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 0)).to_numpy()
    
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('StableStateMat')[['Hinge Number']].count()
    
    more10percent = kappasStSt['amountStSt']>simBound
    kappasStSt = kappasStSt[more10percent]
   
    return kappasStSt

def makeSelectionPureMat(ThisDataPure, ThisFlags, typePureMat):
    
    flagsName = ['NonConv', 'AreaConst']
    possibleFlags = [-1,-5]
    
    selDataMat = ThisFlags.iloc[0,2:5].to_frame().transpose()
    
    matStSt, numMatStSt = np.unique(typePureMat, return_counts= True)
    for i in np.arange(1,5):
        thisStSt = matStSt == i
        if np.sum(thisStSt) == 0:
            selDataMat['Material %d' %i] = 0
        else:
            selDataMat['Material %d' %i] = numMatStSt[thisStSt]
            
    for i in np.arange(5,9):
        thisStSt = matStSt == i
        if np.sum(thisStSt) == 0:
            selDataMat['Non pure material %d' %(i-4)] = 0
        else:
            selDataMat['Non pure material %d' %(i-4)] = numMatStSt[thisStSt]
            
    thisStSt = matStSt == 0   
    if np.sum(thisStSt) == 0:
        selDataMat['Non defined material'] = 0
    else:
        selDataMat['Non defined material'] = numMatStSt[thisStSt]
            
    # selDataMat['amountPure'] = np.size(ThisDataPure,0)
    
    for name, flag in zip(flagsName, possibleFlags):
        numFlags = ThisFlags.iloc[(ThisFlags['Flags'] == flag).values,1].values
        if np.size(numFlags) == 0:
            selDataMat[name] = 0
        else:
            selDataMat[name] = numFlags
            
    return  selDataMat

def makeSelectionVertMat(matstst, allData, Energy, Angles, Curv, tess):
    
    allSpec = pd.DataFrame()
    
    pos = [0,np.int(np.floor(tess[0]/2)),np.int(np.floor(np.prod(tess)/2))]
    
    if tess[0]%2 ==0:
        types = [0]
    else:
        types = np.arange(3)
    
    for j in types:
        i = pos[j]
        specData = SelectPerVertex(matstst[:,i], allData, Energy[:,i], Angles[:,i*4:(i+1)*4], Curv[:,i])
        specData["VertexType"] = j
        if j == 0:
            specData["Neighbours"] = tess[0]-1
        else:
            specData["Neighbours"] = (tess[0]-1)/2
        
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
    selStSt['StdMatEnergy'] = error['TotalEnergy'].values
    selStSt['StdMatCurv'] = error['Curvature'].values
        
    selStSt['Hinge Number'] =  allData.groupby('StableStates')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 0)).to_numpy()
    
    ####Only select stable states that have more than 10% appearance in each simulation
    selStSt['amountStSt'] = allData.groupby('StableStates')[['Hinge Number']].count()
    # more10percent = kappasStSt['amountStSt']>simBound
    # kappasStSt = kappasStSt[more10percent]
    
    
    return selStSt

def deleteNonClustered(selData,ThisFlags):
    
    converg = (selData['StableStates'] != -1).values
    
    clusData = selData.iloc[converg,:]
    
    if (~converg).any():
        non_clustered = selData.iloc[~converg,[-2,-1,0,1]].values.flatten()
        non_clustered[0] = -4
        
        ThisFlags.iloc[(ThisFlags['Flags']==1).values,1] -= non_clustered[1]
        ThisFlags.loc[len(ThisFlags)]=non_clustered
    
    return clusData, ThisFlags

def deleteNonClustered2(allDesigns, allFlags):
    
    mask = (allDesigns['StableStateAll'] != -1)
    
    allDesigns = allDesigns.round(8)
    allFlags = allFlags.round(8)
    kappas = np.unique(allDesigns['kappa'])
    thetas = np.unique(allDesigns['restang'])
    
    for k in kappas:
        for t in thetas:
            thisDesBool = (allDesigns['kappa']==k) & (allDesigns['restang']==t)
            thisConvFlag = (allFlags['kappa']==k) & (allFlags['restang']==t) & (allFlags['Flags'] == 1)
            
            thisDes = allDesigns[thisDesBool]
            
            herenonconv = (thisDes['StableStateAll'] == -1).values
            
            if (herenonconv).any():
                non_clustered = thisDes.iloc[herenonconv,[-1,-2,0,1]].values.flatten()
                non_clustered[0] = -4
                
                if np.sum(thisConvFlag) != 1:
                    print('Theres a problem with the flag counting\n')
                allFlags.iloc[thisConvFlag.values,1] -= non_clustered[1]
                allFlags.loc[len(allFlags)]=non_clustered
    
    return mask, allFlags

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

def SaveForPlotMat(allMat, folder):
    
    possTes = np.unique(allMat['tes'])

    for i in possTes:
        thisfolder = folder + '%d_%d_' %(i,i)
        
        thisMatBool = allMat.iloc[:,2] == i
        thisMat = allMat[thisMatBool] 

        SaveForPlot(thisMat, thisfolder)
    
    return
