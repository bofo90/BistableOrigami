# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:04:51 2019

@author: iniguez
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matl
from mpl_toolkits import mplot3d
import configparser
import os.path
import copy
matl.rcParams['pdf.fonttype'] = 42
matl.rcParams['ps.fonttype'] = 42
matl.rcParams['font.family'] = 'sans-serif'
matl.rcParams['font.sans-serif'] = 'Arial'
matl.rcParams['mathtext.fontset'] = 'cm'


def cm2inch(value):
    return value/2.54

def orderAngles(angles, ucsize, simulations):
    
    symetries = ["" for x in range(simulations)]
    finalAngles = np.zeros((simulations,np.size(angles,1)))
    for sim in np.arange(simulations):
        sym = ''
        ver = np.array(angles[sim,:])
        if np.sum(np.sign(np.around(ver, decimals = 4))) < 0: ### sort to have mayority positive angles
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
        
    return [finalAngles, symetries]
    # return [angles, symetries]

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
    
    matUCstandStSt = np.zeros((np.size(allUCang,0),4))
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

def NiceGraph2D(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
                buffer = [0.0, 0.0, 0.0]):
    gray = '0.2'
    matl.rcParams.update({'font.size': 9})

    if ~np.isnan(mincoord[0]) and ~np.isnan(maxcoord[0]):
        axes.set_xlim([mincoord[0]-buffer[0], maxcoord[0]+buffer[0]])
        if isinstance(divisions[0], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[0]).any():
                axes.set_xticks(divisions[0])
        else:
            if ~np.isnan(divisions[0]):
                axes.set_xticks(np.linspace(mincoord[0],maxcoord[0],divisions[0]))
    axes.set_xlabel(nameX,labelpad=-3, color = gray)
    
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
    axes.set_ylabel(nameY,labelpad=-3, color = gray)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x', colors=gray, direction = 'in', width = 0.4)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray, direction = 'in', width = 0.4)
    axes.tick_params(pad = 2)
    
    axes.tick_params(axis='y', which='minor', colors=gray, direction = 'in', width = 0.4)
    axes.tick_params(axis='x', which='minor', colors=gray, direction = 'in', width = 0.4)
    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
        axes.spines[axis].set_color(gray)
        
    return

def NiceGraph2Dlog(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
                buffer = [0.0, 0.0, 0.0]):
    gray = '0.2'
    matl.rcParams.update({'font.size': 9})

    if ~np.isnan(mincoord[0]) and ~np.isnan(maxcoord[0]):
        axes.set_xlim([mincoord[0]*((buffer[0])**(-1)), maxcoord[0]*(buffer[0])])
        if isinstance(divisions[0], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[0]).any():
                axes.set_xticks(divisions[0])
        else:
            if ~np.isnan(divisions[0]):
                axes.set_xticks(np.logspace(np.log10(mincoord[0]),np.log10(maxcoord[0]),divisions[0]))
#    axes.set_xscale('log')
            
    axes.set_xlabel(nameX)
    if ~np.isnan(mincoord[1]) and ~np.isnan(maxcoord[1]):
        axes.set_ylim([mincoord[1]-buffer[1], maxcoord[1]+buffer[1]])
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.linspace(mincoord[1],maxcoord[1],divisions[1]))
                            
    axes.set_ylabel(nameY)
   
    axes.xaxis.label.set_color(gray)
    axes.tick_params(axis='x',colors=gray)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    axes.tick_params(pad = 2)
    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
        
    return

#%%
def ReadMetadata(file):
    
    if os.path.isfile(file):
        metadata = configparser.RawConfigParser(allow_no_value=True)
        metadata.read(file)
    
        restang = float(metadata.get('options','restang'))
        designang = np.array(metadata.get('options','angDesign').split(),dtype = float)
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

    return restang, designang

plt.close('all')

kappas = np.logspace(-3,1,81)
#kappas = np.logspace(-3,1,13)#23

Folder_name = "Results/Tessellation4/25/2CFF/09-Jan-2020_5_5_"
file_name1 = "/EnergyData.csv" 
file_name2 = "/Hinges.csv"
file_name3 = "/PosStad.csv"
file_name4 = "/Angles.csv"

allDesigns = pd.DataFrame()
allKappasAnalysis = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    
tessellation = np.array(Folder_name.split('_')[1:-1]).astype(int)
numUC = np.prod(tessellation)
  
for subdir in os.listdir(Folder_name):
    if subdir == 'Images':
        continue    
    
    allData = pd.DataFrame()
    
    plt.close('all')
    
    for subdir2 in os.listdir(Folder_name+'/'+subdir):

        folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
    
        dataEnergy = pd.read_csv(folder_name+file_name1)
        simLen = np.size(dataEnergy['Hinge Number'])
        
        dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
        dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']
        dataEnergy['kappa'] = np.ones(simLen)*np.float(subdir2.split('_')[1])  ### '/4' due to error in the computation of the energy
        
        dataCurv = np.loadtxt(folder_name+file_name2,skiprows=1, delimiter = ',', dtype = np.float64)
        dataCurv = np.delete(dataCurv, 0, 1)
        # dataEnergy['Curvature'] = dataHinges['Curvature']
        
        dataLen = np.loadtxt(folder_name+file_name3,skiprows=1, delimiter = ',', dtype = np.float64)
        dataLen = np.delete(dataLen, 0, 1)
        # dataPos = dataPos.reset_index()
        # dataPos = dataPos.rename(columns={'level_0':'Hinge Number', 'level_1':'face1', 'level_2':'face2', 'Hinge Number':'face3', 'Face Areas':'face4'})
        # dataEnergy[['face1','face2','face3','face4']] = dataPos[['face1','face2','face3','face4']]
    
        dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
        dataAngles = np.delete(dataAngles, 0, 1)
        # [dataAnglesOrd, dataSym] = orderAngles(dataAngles, 4, simLen)
        dataAnglesOrd = orderAnglesMat(dataAngles, 4, tessellation)
        # dataStSt = countStableStates(np.resize(dataAnglesOrd,(simLen*numUC,4)), 0.8, 'centroid',True)
        dataEnergy['StableStates'] = countStableStates(dataAnglesOrd, 0.6, 'centroid')
        # dataStSt = np.resize(dataStSt,(simLen,numUC))
                
        allAngles = np.concatenate((dataCurv,dataLen,dataAnglesOrd),1)
        allAngles = pd.DataFrame(allAngles)
        dataEnergy = pd.concat([dataEnergy, allAngles], axis=1, sort=False)
        allData = allData.append(dataEnergy)
        
    restang = np.float(subdir.split('_')[1])
         
    TotSimul = allData.shape[0]
    
    
    minFaceFlag = (allData.iloc[:,10+numUC:-(4*numUC)]<0.11).any(axis = 1)
    allData['Flags'][minFaceFlag] = -5
    
    print(allData['Flags'].value_counts())
    exfl = allData.Flags.values.reshape(TotSimul,1)
    flagmask = (exfl !=1) & (exfl !=2)
    flagmask = ~flagmask.any(axis= 1)
    allData['Mask'] = flagmask
    allData = allData.iloc[flagmask,:]
    
    allData.reset_index(level=0, drop=True)
    
    # [allData['StableStates'], angord] = getStableStatesMat(allData.iloc[:,-1-numUC:-1].to_numpy(), tessellation)
    # orderedang = np.zeros([TotSimul,numUC*4])
    # for i in np.arange(TotSimul):
    #     for j in np.arange(numUC):
    #         orderedang[i,j*4:j*4+4] = allData.iloc[i,-(2+5*numUC)+angord[i,j]*4:-(2+5*numUC)+angord[i,j]*4+4]
    # allData.iloc[:,-(2+5*numUC):-(2+numUC)] = orderedang
            
            
    allData['Curvature'] = allData.iloc[:,10:10+numUC].mean(axis = 1)
    allData['minFace'] = allData.iloc[:,10+numUC:-(2+4*numUC)].min(axis= 1)
    
    colloc = np.concatenate(([1,4,7],np.arange(-(3+4*numUC),-3),[-2,-1]))
    kappasStSt = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates').apply(lambda _df: _df.iloc[:,colloc].mean()))
    
    # kappasStSt['amountSym'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['Symmetry']].nunique())
    # kappasStSt['Symetries'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates').apply(lambda x: str(np.sort(x.Symmetry.unique()))))
    
    kappasStSt['Hinge Number'] =  allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['Hinge Number']].apply(lambda _df2: _df2.sample(1, random_state = 2))).to_numpy()
  
    kappasStSt['restang'] = np.ones(np.shape(kappasStSt)[0])*restang
    
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['DiagonalEnergy']].count())
    more10percent = kappasStSt['amountStSt']>simLen*0.05
    kappasStSt = kappasStSt[more10percent]
    
    kappasStSt = kappasStSt.reset_index(level=0)
    kappasStSt = kappasStSt.reset_index(level=0)#, drop=True)
#    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3','ang4']],1, 'centroid')
#    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3']],0.45, 'centroid')

    kappasnumStSt = kappasStSt.groupby('kappa')['StableStates'].nunique()
    kappasnumStSt = kappasnumStSt.reset_index()
    kappasnumStSt['restang'] = np.ones(np.shape(kappasnumStSt)[0])*restang

    allKappasAnalysis = allKappasAnalysis.append(kappasnumStSt)
    allDesigns = allDesigns.append(kappasStSt)    
      
allDesigns = allDesigns.reset_index(level=0, drop =True)
restangles = allDesigns.restang.drop_duplicates().values


#%%
posStSt = np.shape(allDesigns)[0]
allUCang = np.resize(allDesigns.iloc[:,5:5+numUC*4].values,(posStSt*numUC,4))
[allUCang,sym] = orderAngles(allUCang, 4, posStSt*numUC)
matUCStSt = countStableStates(allUCang, 0.7, 'centroid')
matUCStSt = standarizeStableStates(matUCStSt, allUCang, onlysign = True)
matUCStSt = np.reshape(matUCStSt, (posStSt, numUC))

allDesigns['StableStateAll'] = countStableStates(allDesigns.iloc[:,5:5+numUC*4].values, 1, 'centroid',True)

stst = np.unique(allDesigns['StableStateAll'])
ststUC = np.unique(matUCStSt)

#### Mask results out from the specified range of kappas
# kappasrange = [10**-3,10**0]
# dropers = ~((allDesigns.kappa < kappasrange[0]) | (allDesigns.kappa > kappasrange[-1]))
# allDesigns = allDesigns[dropers]

# cmap2 = matl.cm.get_cmap('Set2',np.size(stst))
cmap2 = matl.cm.get_cmap('jet',np.size(stst))
colors = cmap2(np.linspace(0,1,np.size(stst)))

# cmapUC = matl.cm.get_cmap('Set2',np.size(ststUC))
cmapUC = matl.cm.get_cmap('jet',np.size(ststUC))
colorsUC = cmapUC(np.linspace(0,1,np.size(ststUC)))


#%%

# symstst = ["" for x in stst]

# for i in np.arange(np.size(stst,0)):
#     allsym = allDesigns[allDesigns['StableStateAll'] == stst[i]]['Symetries'].unique()
#     redsym = []
#     for sym in allsym:
#         redsym = np.append(redsym, sym[1:-1].split())
#     symstst[i] = ' '.join(np.unique(redsym))
#     print('Stable state',stst[i],'has symmetries:', symstst[i])
        
#%%
# plt.close('all')   

for i in np.arange(np.size(restangles)):
    
    fig1 = plt.figure(figsize=(cm2inch(4.3), cm2inch(3.1)))
    ax1 = plt.subplot(111)
    fig1.subplots_adjust(top=0.982,
    bottom=0.23,
    left=0.225,
    right=0.940)
    
    thisangBool = allDesigns['restang'] == restangles[i]
    thisang = allDesigns[thisangBool]    
    
    maxCurv = 0.7
    minCurv = 0
    
    NiceGraph2D(ax1, r'$\kappa$', r'$Area$')
    
    ax1.set_ylim([minCurv, maxCurv])
    
    ax1.set_xscale('log')
    ax1.set_xticks([0.001,0.01,0.1,1])
    ax1.set_xlim([0.0007,1.5])
        
    for k in np.arange(np.size(stst)):
        thisstst = thisang[thisang.StableStateAll == stst[k]]
        
        ax1.scatter(thisstst['kappa'], thisstst['minFace'], c = [colors[k-1]], s = 2)

    ax1.axhline(y=0.1, color='r', linestyle='-', linewidth = '0.4')
    
    fig1.show()
    fig1.savefig(Folder_name + '/Images/AreaFaces_' + str(restangles[i].astype(float))+'.pdf', transparent = True)
    fig1.savefig(Folder_name + '/Images/AreaFaces_' + str(restangles[i].astype(float))+ '.png', transparent = True)

#%%
plt.close('all')   

markers = np.array(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])

# yticks = [[-1.2,0,1.0],[-38,-21,-9,0,10],[-2.4,-0.5,0]]
# ylim = [[-1.5,1.5],[-40,13],[-2.7,0.2]]


for i in np.arange(np.size(restangles)):
    
    fig1 = plt.figure(figsize=(cm2inch(4.7), cm2inch(3.3)))
    ax1 = plt.subplot(111)
    fig1.subplots_adjust(top=0.995,
    bottom=0.22,
    left=0.20,
    right=0.98)
    
    thisangBool = allDesigns['restang'] == restangles[i]
    thisang = allDesigns[thisangBool]    
    
    maxTotEn = thisang['TotalEnergy'].max()
    maxCurv = np.ceil(thisang['Curvature'].max())#quantile(0.95)
    minCurv = np.floor(thisang['Curvature'].min())
    
#    NiceGraph2D(ax1, r'$\kappa$', r'$E_\mathrm{tot}$')
    NiceGraph2D(ax1, r'$\kappa$', r'$K$')#, mincoord=[np.log10(kappas[0]),-3], maxcoord=[np.log10(kappas[-1]),3], divisions=[5, 5], buffer=[0.1, 1])
    
# #    ax1.set_yticks(np.linspace(0, maxTotEn,5))
# #    ax1.set_ylim([-0.005, maxTotEn+0.005])
# #    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2f'))
#     ax1.set_yticks(yticks[i])
#     ax1.set_ylim(ylim[i])
# #    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.1f'))
    ax1.set_xscale('log')
#     ax1.set_xticks([0.001,0.01,0.1,1])
#     ax1.set_xlim([0.0007,1.5])
        
    for j in stst:
        thisstst = thisang[thisang.StableStateAll == j]
#    plt.close('all')

#        ax1.scatter(thisstst['kappa'], thisstst['TotalEnergy'], c = colors[thisstst['StableStateAll'].values.astype('int')-1], s = 10)
        ax1.scatter(thisstst['kappa'], thisstst['Curvature'], c = colors[thisstst['StableStateAll']-1], s = 1)#,
#                s = 5, marker = markers[i])#, linestyle = lines[i], lw = 2.5)

#leg = ax1.legend(loc = 2, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False) 
##           borderpad = 0.3, labelspacing = 0.1, handlelength = 0.4, handletextpad = 0.4)
#plt.setp(leg.get_texts(), color='0.2')
#leg.get_frame().set_linewidth(0.4)
    
    fig1.show()
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+'.pdf', transparent = True)
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+ '.png', transparent = True)
    fig1.savefig(Folder_name + '/Images/Restang_' + str(restangles[i].astype(float))+'.pdf', transparent = True)
    fig1.savefig(Folder_name + '/Images/Restang_' + str(restangles[i].astype(float))+ '.png', transparent = True)

#%%
fig3 = plt.figure(figsize=(cm2inch(10), cm2inch(7)))
ax3 = plt.subplot(111,projection='3d')
ax3.set_xlim([-np.pi,np.pi])
ax3.set_ylim([-np.pi,np.pi])
ax3.set_zlim([-np.pi,np.pi])

markers = np.array(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
    
# ax3.scatter(allUCang[:,0],allUCang[:,1],allUCang[:,2], c = colorsUC[matUCStSt.flatten()-1])
# for i in ststUC:
#     ax3.scatter(-400,-400,-400, c = [colorsUC[i-1]], label = i)

for i, j in zip(np.arange(numUC), markers):
    ax3.scatter(allUCang[i::numUC,0],allUCang[i::numUC,1],allUCang[i::numUC,2], 
                c = colors[allDesigns['StableStateAll'].values-1], marker = j)
for i in stst:
    ax3.scatter(-400,-400,-400, c = [colors[i-1]], label = i)

plt.legend()

#%%
allDesigns[['kappa','Hinge Number','StableStateAll','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforAllImages.csv', index = False)
allDesigns.groupby('StableStateAll').apply(lambda df: df.sample(1, random_state = 0))[['kappa','Hinge Number','StableStateAll','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforStableStatesImages.csv', index = False)