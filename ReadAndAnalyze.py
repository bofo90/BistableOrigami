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
    
    numhin, numedg, restang, x, y = ReadMetadataMat(folder_name+file_metadata)
    numUC = np.int(x*y)
    
    dataParam = pd.read_csv(folder_name+file_name1)
    simLen = np.size(dataParam['Hinge Number'])   
    
    metadata = pd.DataFrame()
    [ang, kap] = getRAandK(folder_name)
    metadata['kappa'] = np.ones(simLen)*kap
    if oldSample(folder_name):
        metadata['kappa'] = metadata['kappa']/4
    metadata['restang'] = np.ones(simLen)*ang/np.pi
    metadata['tes'] = np.int(y)
        
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
    # Energy = Energy/(kap*4+8)
    
    # Energy = np.sum(LenEnergy, axis = 2)/8+np.sum(AngEnergy, axis = 2)/(4*kap)
    
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
    dataEnergy['TotalEnergy'] = dataEnergy['TotalEnergy']# /(kap*numhin+numedg)
    dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']/kap/numhin
    dataEnergy['EdgeEnergy'] = dataEnergy['EdgeEnergy']/numedg
    # dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['HingeEnergy']
    
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
    
    #### assign flag 1 to converged simulations and 
    ThisData.loc[~flagmask, 'Flags'] = 1
    
    #### assign -5 to converged results with min area
    onlyminArea = minFaceFlag & ~flagmask  
    ThisData.loc[onlyminArea,'Flags'] = -5
    
    #### assign -1 to non-converged ones
    ThisData.loc[flagmask, 'Flags'] = -1
    
    
    flagStad = ThisData['Flags'].value_counts()
    flagStad = flagStad.reset_index(level=0)
    flagStad = flagStad.rename(columns={"index": "Flags", "Flags": "amountFlags"})
    flagStad['kappa'] = np.ones((np.size(flagStad,0),1))*ThisData.iloc[0,0]
    flagStad['restang'] = np.ones((np.size(flagStad,0),1))*ThisData.iloc[0,1]
    
    #add the flag -5 to the flag 1 since it is considered converged
    here_minFace = flagStad.iloc[:,0].values == -5
    here_conv = flagStad.iloc[:,0].values == 1
    if here_minFace.any() & here_conv.any():
        flagStad.iloc[here_conv,1] += flagStad.iloc[here_minFace,1].values
    if here_minFace.any() & ~here_conv.any():
        flagStad = flagStad.append(flagStad.iloc[here_minFace,:], ignore_index = True)
        flagStad.iloc[-1,0] = 1
    
    #### mask the results 
    convmask = ((ThisData.Flags.values) == 1) | ((ThisData.Flags.values) == -5)
    MaskedData = ThisData.iloc[convmask,:]
    MaskedData = MaskedData.drop(['Flags'], axis = 1)
    # MaskedData = ThisData
    
    if returnMask & ~returnStad:
        return MaskedData, convmask
    if returnStad & ~returnMask:
        return MaskedData, flagStad
    if returnStad & returnMask:
        return MaskedData, convmask, flagStad
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
    
    for i in np.arange(10,14):
        thisstate = inverse == i
        inverse[thisstate] = np.argmin(distances[thisstate,:10],axis = 1)
    
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
    matType = np.zeros((simulations,3))
    
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
        
        #check for vert lines pattern with miura
        convMat = signal.convolve2d(simStSt[i,:,:], convMat2, mode = 'valid')
        matType[i,1] = np.sum((convMat == 0) & (sumConMat == 18))/numUCmat
        
        #check for horz lines pattern with miura
        convMat = signal.convolve2d(simStSt[i,:,:], convMat3, mode = 'valid')
        matType[i,1] += np.sum((convMat == 0) & (sumConMat == 10))/numUCmat

        #check for all same with fold from rigid origami
        matType[i,2] = np.sum(sumConMat == 32)/numUCmat
        matType[i,2] += np.sum(sumConMat == 36)/numUCmat
   
    purity = np.max(matType, axis = 1) == 1 
   
    matError = np.zeros((simulations,3))
    mat1Error = np.zeros((simulations,1))
    mat2Error = np.zeros((simulations,6))
    
    convMat5 = np.array([[1,2],[-1,2]])
    convMat6 = np.array([[1,10],[1,0]]) 
    convMat7 = np.array([[1,1],[10,0]])
    convMat8 = np.array([[0,0],[0,1]])
    convMat9 = np.array([[1,-1],[0,0]])
    convMat10 = np.array([[1,0],[-1,0]])
    
    
    stableStateVert = np.zeros((simulations, 15))
    
    for i in np.arange(simulations):
        
        stver, numstst = np.unique(simStSt[i,:,:], return_counts = True)
        stableStateVert[i,stver] = numstst
        
        if purity[i]:
            continue

        sumConMat = signal.convolve2d(simStSt[i,:,:], convMat4, mode = 'valid')
        
        #### Material 1 defects
        #check for dislocations in mat 1
        convMat = signal.convolve2d(simStSt[i,:,:], convMat2, mode = 'valid')
        mat1Error[i] += np.sum((convMat == 0) & (sumConMat == 2))
        matError[i, 0] += np.sum((convMat == 0) & (sumConMat == 2))/numUCmat
        convMat = signal.convolve2d(simStSt[i,:,:], convMat3, mode = 'valid')
        mat1Error[i] += np.sum((convMat == 0) & (sumConMat == 2))
        matError[i, 0] += np.sum((convMat == 0) & (sumConMat == 2))/numUCmat
        
        #check for crossing of dislocations in mat 1
        mat1Error[i] += np.sum((convMat == 0) & (sumConMat == 0))*2
        matError[i,0] += np.sum((convMat == 0) & (sumConMat == 0))/numUCmat
        mat1Error[i] += np.sum((convMat == 0) & (sumConMat == 4))*2
        matError[i,0] += np.sum((convMat == 0) & (sumConMat == 4))/numUCmat
        
        #check for flat state
        matError[i,2] += np.sum((convMat == 0) & (sumConMat == 56))/numUCmat
        
        #### Material 2 defects
        #check for dislocation in mat 2
        convMat = signal.convolve2d(simStSt[i,:,:], convMat1, mode = 'valid')
        matError[i,1] += np.sum((convMat == 0) & (sumConMat == 10))/numUCmat
        mat2Error[i,2] += np.sum((convMat == 0) & (sumConMat == 10))
        matError[i,1] += np.sum((convMat == 0) & (sumConMat == 18))/numUCmat
        mat2Error[i,2] += np.sum((convMat == 0) & (sumConMat == 18))
        
        #check for other type of dislocation in mat 2
        matError[i,1] += np.sum((convMat == 0) & (sumConMat == 8))/numUCmat
        matError[i,1] += np.sum((convMat == 0) & (sumConMat == 12))/numUCmat
        matError[i,1] += np.sum((convMat == 0) & (sumConMat == 16))/numUCmat
        matError[i,1] += np.sum((convMat == 0) & (sumConMat == 20))/numUCmat
        mat2Error[i,3] += np.sum((convMat == 0) & (sumConMat == 8))
        mat2Error[i,3] += np.sum((convMat == 0) & (sumConMat == 12))
        mat2Error[i,3] += np.sum((convMat == 0) & (sumConMat == 16))
        mat2Error[i,3] += np.sum((convMat == 0) & (sumConMat == 20))
        
        #check for change from mat 2 to mat 3 or to mat 2 to another direction
        convMat = signal.convolve2d(simStSt[i,:,:], convMat5, mode = 'valid')
        matError[i,1] += np.sum((convMat == 10) & (sumConMat == 23))/numUCmat/2
        matError[i,2] += np.sum((convMat == 10) & (sumConMat == 23))/numUCmat/2
        mat2Error[i,0] += np.sum((convMat == 10) & (sumConMat == 23))
        
        matError[i,1] += np.sum((convMat == 10) & (sumConMat == 15))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 10) & (sumConMat == 15))
        matError[i,1] += np.sum((convMat == 10) & (sumConMat == 13))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 10) & (sumConMat == 13))
        
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat5,2), mode = 'valid')
        matError[i,1] += np.sum((convMat == 10) & (sumConMat == 23))/numUCmat/2
        matError[i,2] += np.sum((convMat == 10) & (sumConMat == 23))/numUCmat/2
        mat2Error[i,0] += np.sum((convMat == 10) & (sumConMat == 23))
        
        matError[i,1] += np.sum((convMat == 10) & (sumConMat == 15))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 10) & (sumConMat == 15))
        matError[i,1] += np.sum((convMat == 10) & (sumConMat == 13))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 10) & (sumConMat == 13))
        
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat5,1), mode = 'valid')
        matError[i,1] += np.sum((convMat == 18) & (sumConMat == 25))/numUCmat/2
        matError[i,2] += np.sum((convMat == 18) & (sumConMat == 25))/numUCmat/2
        mat2Error[i,0] += np.sum((convMat == 18) & (sumConMat == 25))
        
        matError[i,1] += np.sum((convMat == 18) & (sumConMat == 15))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 18) & (sumConMat == 15))
        matError[i,1] += np.sum((convMat == 18) & (sumConMat == 13))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 18) & (sumConMat == 13))
        
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat5,3), mode = 'valid')
        matError[i,1] += np.sum((convMat == 18) & (sumConMat == 25))/numUCmat/2
        matError[i,2] += np.sum((convMat == 18) & (sumConMat == 25))/numUCmat/2
        mat2Error[i,0] += np.sum((convMat == 18) & (sumConMat == 25))
        
        matError[i,1] += np.sum((convMat == 18) & (sumConMat == 15))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 18) & (sumConMat == 15))
        matError[i,1] += np.sum((convMat == 18) & (sumConMat == 13))/numUCmat
        mat2Error[i,4] += np.sum((convMat == 18) & (sumConMat == 13))
        
        #check for corners between two types of miura-ori
        convMatCorner = signal.convolve2d(simStSt[i,:,:], convMat8, mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], convMat6, mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))
        convMat = signal.convolve2d(simStSt[i,:,:], convMat7, mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))
        
        convMatCorner = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,2), mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat6,2), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat7,2), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))
        
        convMatCorner = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,1), mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat6,1), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat7,1), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))
        
        convMatCorner = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,3), mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat6,3), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 26) | (convMat == 34)))
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat7,3), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 58) | (convMat == 50)))
        
        #check for corners with dislocations
        convMatCorner = signal.convolve2d(simStSt[i,:,:], convMat8, mode = 'valid')
        convMatEqual = signal.convolve2d(simStSt[i,:,:], convMat9, mode = 'valid')
        convMatEqual2 = signal.convolve2d(simStSt[i,:,:], convMat10, mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,2), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        
        convMatCorner = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,2), mode = 'valid')
        convMatEqual = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat9,2), mode = 'valid')
        convMatEqual2 = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat10,2), mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,3), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        
        convMatCorner = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,1), mode = 'valid')
        convMatEqual = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat9,1), mode = 'valid')
        convMatEqual2 = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat10,1), mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,2), mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        
        convMatCorner = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,3), mode = 'valid')
        convMatEqual = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat9,3), mode = 'valid')
        convMatEqual2 = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat10,3), mode = 'valid')
        convMat = signal.convolve2d(simStSt[i,:,:], convMat8, mode = 'valid')
        matError[i,1] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 2) | (convMatCorner == 3)) & ((convMat == 4) | (convMat == 5)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        matError[i,1] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))/numUCmat
        mat2Error[i,5] += np.sum(((convMatCorner == 4) | (convMatCorner == 5)) & ((convMat == 2) | (convMat == 3)) & (convMatEqual == 0) & (convMatEqual2 == 0))
        
        #check for change from mat 2 to mat 3
        convMatCorner = signal.convolve2d(simStSt[i,:,:], convMat8, mode = 'valid')
        convMatCorner2 = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat8,2), mode = 'valid')
        
        convMatEqual = signal.convolve2d(simStSt[i,:,:], convMat9, mode = 'valid')
        convMatEqual2 = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat9,2), mode = 'valid')
        pos = np.sum((convMatEqual == 0) & (convMatEqual2 == 0) & 
                    ((convMatCorner == 2) | (convMatCorner == 3) | (convMatCorner == 4) | (convMatCorner == 5)) & 
                    ((convMatCorner2 == 8) | (convMatCorner2 == 9)))
        matError[i,1] += pos/numUCmat/2
        matError[i,2] += pos/numUCmat/2
        mat2Error[i,1] += pos
        pos = np.sum((convMatEqual == 0) & (convMatEqual2 == 0) & 
                    ((convMatCorner2 == 2) | (convMatCorner2 == 3) | (convMatCorner2 == 4) | (convMatCorner2 == 5)) & 
                    ((convMatCorner == 8) | (convMatCorner == 9)))
        matError[i,1] += pos/numUCmat/2
        matError[i,2] += pos/numUCmat/2
        mat2Error[i,1] += pos
        
        convMatEqual = signal.convolve2d(simStSt[i,:,:], convMat10, mode = 'valid')
        convMatEqual2 = signal.convolve2d(simStSt[i,:,:], np.rot90(convMat10,2), mode = 'valid')
        pos = np.sum((convMatEqual == 0) & (convMatEqual2 == 0) & 
                    ((convMatCorner == 2) | (convMatCorner == 3) | (convMatCorner == 4) | (convMatCorner == 5)) & 
                    ((convMatCorner2 == 8) | (convMatCorner2 == 9)))
        matError[i,1] += pos/numUCmat/2
        matError[i,2] += pos/numUCmat/2
        mat2Error[i,1] += pos
        pos = np.sum((convMatEqual == 0) & (convMatEqual2 == 0) & 
                    ((convMatCorner2 == 2) | (convMatCorner2 == 3) | (convMatCorner2 == 4) | (convMatCorner2 == 5)) & 
                    ((convMatCorner == 8) | (convMatCorner == 9)))
        matError[i,1] += pos/numUCmat/2
        matError[i,2] += pos/numUCmat/2
        mat2Error[i,1] += pos
    
    materials = matType + matError
    # if (np.round(np.sum(materials,axis = 1),5) < 0.9).any():
    if (np.round(np.sum(materials,axis = 1),5) != 1).any():
        print('Material detection error')
        here = np.arange(np.size(materials,0))
        here = here[np.round(np.sum(materials,axis = 1),5) != 1]
        
    if (np.round(np.sum(materials,axis = 1),5) > 1).any():
        print('Over assigning material type')
        
    matName = np.argmax(materials, axis = 1)+1
    matName[~purity] += 3
    
    nomat = np.max(materials, axis = 1) == 0
    matName[nomat] = 0
    
    defects = mat1Error / (tessellation[0]-1)
    
    mat1nopure = (np.sum(stableStateVert[:,:2],axis = 1) != np.prod(tessellation)) & (matName == 4)
    defects[mat1nopure] = np.array([np.sum(stableStateVert[mat1nopure,2:],axis = 1) + 9]).T#/(tessellation[0]-1)
    
    mat2nopure = matName == 5
    defects[mat2nopure,0] = (mat2Error[mat2nopure,5]+mat2Error[mat2nopure,4])

    lineError = (1-np.max(matType,1))*numUCmat/(tessellation[0]-1)
    
    return purity, matName, np.max(matType,1), defects

def getPureMatComp(simStSt, tessellation):
    
    numUCmat = np.prod(tessellation-1)
    simulations = np.size(simStSt,0)
    simStSt = simStSt.reshape((simulations,*tessellation))
    simMatAna = np.zeros((simulations, *tessellation-1))
    
    matType = np.zeros((simulations,12))
        
    stableStateVert = np.zeros((simulations, 15))
    
    for i in np.arange(simulations):
        
        stver, numstst = np.unique(simStSt[i,:,:], return_counts = True)
        stableStateVert[i,stver] = numstst
        
        unitcellmat = np.zeros((*tessellation-1,4))
        unitcellshift = np.zeros(tessellation-1)
        for x in np.arange(np.size(simStSt[i,:,:],0)-1):
            for y in np.arange(np.size(simStSt[i,:,:],1)-1):
                unitcell = simStSt[i,x:x+2,y:y+2]
                orunitcell, unitcellshift[x,y] = orderUnitCell(unitcell)
                unitcellmat[x,y,:] = orunitcell.flatten()
        unitcellComb = np.concatenate([unitcellmat.reshape(numUCmat,4), unitcellshift.reshape(numUCmat,1)], axis = 1)
        
        unitcellCombUn, unitcellLoc = np.unique(unitcellComb, axis = 0, return_inverse = True)
        unitcellUn = unitcellCombUn[:,:4]
        unitcellshiftUn = unitcellCombUn[:,4]
        
        materialStableStates = np.zeros(np.size(unitcellUn,0))
        
        #check for main materials
        materialStableStates[(unitcellUn == [0,1,1,0]).all(axis = 1)] = 1
        materialStableStates[(unitcellUn == [2,2,3,3]).all(axis = 1) & (unitcellshiftUn%2 == 0)] = 2
        materialStableStates[(unitcellUn == [4,4,5,5]).all(axis = 1) & (unitcellshiftUn%2 == 1)] = 2
        materialStableStates[(unitcellUn == [8,8,8,8]).all(axis = 1)] = 3
        materialStableStates[(unitcellUn == [9,9,9,9]).all(axis = 1)] = 3
        
        #check for grain boundaries mat1
        #twin
        materialStableStates[(unitcellUn == [0,0,1,1]).all(axis = 1)] = 4
        materialStableStates[(unitcellUn == [1,1,1,1]).all(axis = 1)] = 4
        materialStableStates[(unitcellUn == [0,0,0,0]).all(axis = 1)] = 4
        
        #check for precipitation of mat 2 in mat1 at grain boundaries
        materialStableStates[(unitcellUn[:,0:2]<2).all(axis = 1) & 
                             ((unitcellUn[:,2:]>1).all(axis = 1) & (unitcellUn[:,2:]<6).all(axis = 1))] = 5
        #precipitation at corners
        materialStableStates[(unitcellUn[:,:3]<2).all(axis = 1) &
                             ((unitcellUn[:,3]>1) & (unitcellUn[:,3]<6))] = 5
        materialStableStates[(unitcellUn[:,[0,1,3]]<2).all(axis = 1) &
                             ((unitcellUn[:,2]>1) & (unitcellUn[:,2]<6))] = 5  
        materialStableStates[(unitcellUn[:,0]<2) &
                             ((unitcellUn[:,1:]>1).all(axis = 1) & (unitcellUn[:,1:]<6).all(axis = 1))] = 5
        #precipitation at double corner
        materialStableStates[(unitcellUn[:,[0,3]]<2).all(axis = 1) & 
                             ((unitcellUn[:,1:3]>1).all(axis = 1) & (unitcellUn[:,1:3]<6).all(axis = 1))] = 5
                
        #check for precipitation of mat 3 in mat1 at grain boundaries
        materialStableStates[((unitcellUn[:,:2]<2).all(axis = 1)) &
                             ((unitcellUn[:,2:]>7).all(axis = 1) & (unitcellUn[:,2:]<10).all(axis = 1))] = 6 
        #precipitation at corners
        materialStableStates[((unitcellUn[:,:3]<2).all(axis = 1)) &
                             ((unitcellUn[:,3]>7) & (unitcellUn[:,3]<10))] = 6 
        materialStableStates[((unitcellUn[:,[0,1,3]]<2).all(axis = 1)) &
                             ((unitcellUn[:,2]>7) & (unitcellUn[:,2]<10))] = 6  
        materialStableStates[((unitcellUn[:,0]<2)) &
                             ((unitcellUn[:,1:]>7).all(axis = 1) & (unitcellUn[:,1:]<10).all(axis = 1))] = 6
        
        
        #check for grain boundary mat2
        #twin
        materialStableStates[(unitcellUn == [2,2,2,2]).all(axis = 1)] = 7
        materialStableStates[(unitcellUn == [3,3,3,3]).all(axis = 1)] = 7
        materialStableStates[(unitcellUn == [4,4,4,4]).all(axis = 1)] = 7
        materialStableStates[(unitcellUn == [5,5,5,5]).all(axis = 1)] = 7
        
        #slip
        materialStableStates[(unitcellUn == [2,3,3,2]).all(axis = 1)] = 8
        materialStableStates[(unitcellUn == [4,5,5,4]).all(axis = 1)] = 8
        #corner slip and twin
        materialStableStates[(unitcellUn == [2,2,3,3]).all(axis = 1) & (unitcellshiftUn%2 == 1)] = 8
        materialStableStates[(unitcellUn == [4,4,5,5]).all(axis = 1) & (unitcellshiftUn%2 == 0)] = 8   
        
        #rotation
        materialStableStates[((unitcellUn[:,:2]>1).all(axis = 1) & (unitcellUn[:,:2]<4).all(axis = 1)) &
                             ((unitcellUn[:,2:]>3).all(axis = 1) & (unitcellUn[:,2:]<6).all(axis = 1))] = 9   
        #rotation at corners
        materialStableStates[((unitcellUn[:,:3]>1).all(axis = 1) & (unitcellUn[:,:3]<4).all(axis = 1)) &
                             ((unitcellUn[:,3]>3) & (unitcellUn[:,3]<6))] = 9  
        materialStableStates[((unitcellUn[:,[0,1,3]]>1).all(axis = 1) & (unitcellUn[:,[0,1,3]]<4).all(axis = 1)) &
                             ((unitcellUn[:,2]>3) & (unitcellUn[:,2]<6))] = 9  
        materialStableStates[((unitcellUn[:,0]>1) & (unitcellUn[:,0]<4)) &
                             ((unitcellUn[:,1:]>3).all(axis = 1) & (unitcellUn[:,1:]<6).all(axis = 1))] = 9
        #rotation at double corners
        materialStableStates[((unitcellUn[:,[0,3]]>1).all(axis = 1) & (unitcellUn[:,[0,3]]<4).all(axis = 1)) &
                             ((unitcellUn[:,1:3]>3).all(axis = 1) & (unitcellUn[:,1:3]<6).all(axis = 1))] = 9  
 
        
        #check interface material 2 and 3
        materialStableStates[((unitcellUn[:,:2]>1).all(axis = 1) & (unitcellUn[:,:2]<6).all(axis = 1)) &
                             ((unitcellUn[:,2:]>7).all(axis = 1) & (unitcellUn[:,2:]<10).all(axis = 1))] = 10
        #interfaces at corners
        materialStableStates[((unitcellUn[:,:3]>1).all(axis = 1) & (unitcellUn[:,:3]<6).all(axis = 1)) &
                             ((unitcellUn[:,3]>7) & (unitcellUn[:,3]<10))] = 10
        materialStableStates[((unitcellUn[:,[0,1,3]]>1).all(axis = 1) & (unitcellUn[:,[0,1,3]]<6).all(axis = 1)) &
                             ((unitcellUn[:,2]>7) & (unitcellUn[:,2]<10))] = 10 
        materialStableStates[((unitcellUn[:,0]>1) & (unitcellUn[:,0]<6)) &
                             ((unitcellUn[:,1:]>7).all(axis = 1) & (unitcellUn[:,1:]<10).all(axis = 1))] = 10
        #interface at double corners
        materialStableStates[((unitcellUn[:,[0,3]]>1).all(axis = 1) & (unitcellUn[:,[0,3]]<6).all(axis = 1)) &
                             ((unitcellUn[:,1:3]>7).all(axis = 1) & (unitcellUn[:,1:3]<10).all(axis = 1))] = 9 
        
        #check interface between material 1, 2 and 3
        materialStableStates[(unitcellUn < 2).any(axis = 1) & 
                             ((unitcellUn > 1) & (unitcellUn < 6)).any(axis = 1)  &
                             ((unitcellUn > 7) & (unitcellUn < 10)).any(axis = 1)] = 11
        
        #check if vertex is in impossible state
        materialStableStates[(unitcellUn == 6).any(axis = 1) | (unitcellUn == 7).any(axis = 1)] = 12
        
        
        simMatAna[i,:,:] = materialStableStates[unitcellLoc].reshape(*tessellation-1)
        
        types, typescounts = np.unique(simMatAna[i,:,:], return_counts = True)
        
        if (types == 0).any():
            typescounts = np.delete(typescounts, np.where(types == 0))
            types = np.delete(types, np.where(types == 0))

        matType[i, types.astype(int)-1] = typescounts/numUCmat
        
    
    if (np.round(np.sum(matType,axis = 1),5) != 1).any():
        print('Material detection error')
        here = np.arange(np.size(matType,0))
        here = here[np.round(np.sum(matType,axis = 1),5) != 1]
        
    material = np.zeros((simulations,3))
    material[:,0] = np.sum(matType[:,[0,3]], axis = 1) + np.sum(matType[:,[4,5]], axis = 1)/2 + matType[:,10]/3
    material[:,1] = np.sum(matType[:,[1,6,7,8]], axis = 1) + np.sum(matType[:,[4,9]], axis = 1)/2 + matType[:,10]/3
    material[:,2] = matType[:,2]+ np.sum(matType[:,[5,9]], axis = 1)/2 + matType[:,10]/3
        
    purity = ~(matType[:,3:]>0).any(axis = 1)
    matName = np.argmax(material, axis = 1)+1
    matName[~purity] += 3
    
    nomat = (np.max(material, axis = 1) == 0) | (matType[:,11] != 0)
    matName[nomat] = 0
    
    grainSizes, grainDis = calculateGrainSize(simMatAna)
    
    return purity, matName, np.max(matType[:,:3],1), matType[:,3:-1], grainSizes

def orderUnitCell(unitcell):
    
    shift = np.argmin(unitcell)
    if shift > 1:
        shift = 5-shift
        
    unitcell = np.rot90(unitcell,shift)
    
    if unitcell[1,0] < unitcell[0,1]:
        unitcell = unitcell.T
        shift += 1
        
    return unitcell, shift

def calculateGrainSize(simMatAna):
    
    from sklearn.cluster import DBSCAN
    
    simulations = np.size(simMatAna,0)
    sizes = np.zeros((simulations,3))
    materials = np.arange(3)
    
    x = np.size(simMatAna,1)
    y = np.size(simMatAna,2)
    xPos, yPos = np.meshgrid(np.arange(x),np.arange(y))
    matUCTemp = np.zeros(np.shape(simMatAna[0,:,:]))
    matGrainDis = np.zeros((simulations,x+1,y+1))
    
    for i in np.arange(simulations):
        for m in materials:
            hereMat = simMatAna[i,:,:] == m+1
            
            if ~hereMat.any():
                continue
            
            points = np.concatenate(([xPos[hereMat].flatten()],[yPos[hereMat].flatten()]), axis = 0).T

            db = DBSCAN(eps = 1, min_samples = 1).fit(points)
            inverse = db.labels_+1
            
            matUCTemp[points[:,1],points[:,0]] = inverse
            matGrainDis[i,:,:] = convertToMaterial(matUCTemp)
            matUCTemp = matUCTemp*0
            
            
            grainsizeTempX = 0
            xTemp = x
            for j in np.arange(x):
                grain, grainsize = np.unique(matGrainDis[i,j,:], return_counts = True)
                if (grain == 0).any():
                    grainsize = np.delete(grainsize, np.where(grain == 0))      
                if (grain == 0).all():
                    xTemp -= 1
                else:
                    grainsizeTempX += np.mean(grainsize)
                
            grainsizeTempY = 0
            yTemp = y
            for k in np.arange(y):
                grain, grainsize = np.unique(matGrainDis[i,:,k], return_counts = True)
                if (grain == 0).any():
                    grainsize = np.delete(grainsize, np.where(grain == 0))
                if (grain == 0).all():
                    yTemp -= 1
                else:
                    grainsizeTempY += np.mean(grainsize)
                    
            sizes[i,m] = (grainsizeTempX/xTemp+grainsizeTempY/yTemp)/2
                    
            # grain, grainsize = np.unique(matGrainDis[i,:,:], return_counts = True)
            # if (grain == 0).any():
            #     grainsize = np.delete(grainsize, np.where(grain == 0))
            #     grain = np.delete(grain, np.where(grain == 0))
            # sizes[i,m] = np.mean(grainsize)
    
    return sizes, matGrainDis

def convertToMaterial(matUCTemp):
    
    x = np.size(matUCTemp,0)+1
    y = np.size(matUCTemp,1)+1
    
    matSV = np.zeros((x,y))
    
    for i in np.arange(x-1):
        for j in np.arange(y-1):
            if matUCTemp[i,j] != 0:
                matSV[i:i+2,j:j+2] = matUCTemp[i,j]
    
    return matSV

def getPureMatLine(simStSt, tessellation):
    
    from scipy import signal
    
    numUCmat = np.prod(tessellation-1)
    simulations = np.size(simStSt,0)
    simStSt = simStSt.reshape((simulations,*tessellation))
    matType = np.zeros((simulations,3))
    
    ##### not completed!
    
    return purity, matName, np.max(matType,1), mat1Error
    
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

def makeSelectionPerStSt2(allData):
    
    
    kappasStSt = allData.groupby('StableStateAll').apply(lambda _df: _df.mean())
        
    kappasStSt['Hinge Number'] =  allData.groupby('StableStateAll')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 0)).to_numpy()
    
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('StableStateAll')[['amountStSt']].sum()
   
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
    for i in np.arange(1,4):
        thisStSt = matStSt == i
        if np.sum(thisStSt) == 0:
            selDataMat['Material %d' %i] = 0
        else:
            selDataMat['Material %d' %i] = numMatStSt[thisStSt]
            
    for i in np.arange(4,7):
        thisStSt = matStSt == i
        if np.sum(thisStSt) == 0:
            selDataMat['Non pure material %d' %(i-3)] = 0
        else:
            selDataMat['Non pure material %d' %(i-3)] = numMatStSt[thisStSt]
            
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
        non_clustered[0] = -3
        
        ThisFlags.iloc[(ThisFlags['Flags']==1).values,1] -= non_clustered[1]
        ThisFlags.loc[len(ThisFlags)]=non_clustered
    
    return clusData, ThisFlags

def makeSelectionClustered(allDesigns, allFlags, simTresh):
    
    allDesignsRed = pd.DataFrame()
    
    allDesigns = allDesigns.round(8)
    allFlags = allFlags.round(8)
    kappas = np.unique(allDesigns['kappa'])
    thetas = np.unique(allDesigns['restang'])
    
    for k in kappas:
        for t in thetas:
            thisDesBool = (allDesigns['kappa']==k) & (allDesigns['restang']==t)
            thisConvFlag = (allFlags['kappa']==k) & (allFlags['restang']==t) & (allFlags['Flags'] == 1)
            
            thisDes = allDesigns[thisDesBool]
            
            thisDesRed = makeSelectionPerStSt2(thisDes)
            
            herenonconv = (thisDesRed['StableStateAll'] == -1).values
            
            if (herenonconv).any():
                non_clustered = thisDesRed.iloc[herenonconv,[-1,-2,0,1]].values.flatten()
                non_clustered[0] = -4
                
                if np.sum(thisConvFlag) != 1:
                    print('Theres a problem with the flag counting\n')
                allFlags.iloc[thisConvFlag.values,1] -= non_clustered[1]
                allFlags.loc[len(allFlags)]=non_clustered
                
            herenonconv = ((thisDesRed['StableStateAll'] != -1) & (thisDesRed['amountStSt']<= simTresh)).values
            
            if (herenonconv).any():
                non_clustered_all = thisDesRed.iloc[herenonconv,[-1,-2,0,1]].values
                non_clustered = np.sum(non_clustered_all,axis = 0)
                non_clustered[0] = -6
                non_clustered[[2,3]] = non_clustered_all[0,[2,3]]
                
                if np.sum(thisConvFlag) != 1:
                    print('Theres a problem with the flag counting\n')
                allFlags.iloc[thisConvFlag.values,1] -= non_clustered[1]
                allFlags.loc[len(allFlags)]=non_clustered
                
            
            allDesignsRed = allDesignsRed.append(thisDesRed)
    
    allDesignsRed = allDesignsRed.reset_index(level=0, drop =True)
    
    mask = (allDesignsRed['StableStateAll'] == -1) | (allDesignsRed['amountStSt'] <= simTresh)
    
    return allFlags, allDesignsRed[~mask]

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

    return hing, edge, restang, x, y

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
    
    des[['kappa','Hinge Number','StableStateMat','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforAllImages.csv', index = False)
    des.groupby('StableStateMat').apply(lambda df: df.sample(1, random_state = 0))[['kappa','Hinge Number','StableStateMat','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforStableStatesImages.csv', index = False)
       
    return

def SaveForPlotMat(allMat, folder):
    
    possTes = np.unique(allMat['tes'])

    for i in possTes:
        thisfolder = folder + '%d_%d_' %(i,i)
        
        thisMatBool = allMat.iloc[:,2] == i
        thisMat = allMat[thisMatBool] 

        SaveForPlot(thisMat, thisfolder)
    
    return
