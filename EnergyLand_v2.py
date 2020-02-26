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

Folder_name = "Results/SingleVertex4/2CFF/02-Dec-2019_0.00_90.00_180.00_270.00_"

allDesigns = pd.DataFrame()
# allKappasAnalysis = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    
designang = np.array(Folder_name.split('_')[1:-1]).astype(float)
numvertex = np.int(Folder_name.split('/')[1][-1])
  
for subdir in os.listdir(Folder_name):
    if subdir == 'Images':
        continue      
    
    for subdir2 in os.listdir(Folder_name+'/'+subdir):
        folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
        
        ThisData, simLen = raa.ReadFile(folder_name)
        ThisData = raa.maskBadResults(ThisData, True) 
        
        simLen = np.size(ThisData,0)
        ThisData.iloc[:,-numvertex::],sym = raa.orderAngles(ThisData.iloc[:,-numvertex::].to_numpy(), numvertex, simLen)
        ThisData['StableStates'] = raa.countStableStates(ThisData.iloc[:,-numvertex::], 0.6, 'centroid')
        
        selData = raa.makeSelectionPerStSt(ThisData, simLen*0.0)
        
        allDesigns = allDesigns.append(selData)
        
    # kappasStSt['amountSym'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['Symmetry']].nunique())
    # kappasStSt['Symetries'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates').apply(lambda x: str(np.sort(x.Symmetry.unique()))))
    
    # kappasnumStSt = kappasStSt.groupby('kappa')['StableStates'].nunique()
    # kappasnumStSt = kappasnumStSt.reset_index()
    # kappasnumStSt['restang'] = np.ones(np.shape(kappasnumStSt)[0])*restang

    # allKappasAnalysis = allKappasAnalysis.append(kappasnumStSt)
    
allDesigns = allDesigns.reset_index(level=0, drop =True)
allDesigns = allDesigns.round(8)
restangles = allDesigns.restang.drop_duplicates().values
kappas = allDesigns.kappa.drop_duplicates().values

#%%
allDesAng = allDesigns.iloc[:,8:8+numvertex].values
allDesigns['StableStateAll'] = raa.countStableStates(allDesAng, 0.9, 'centroid', True)
allDesigns['StableStateAll'] = raa.standarizeStableStates(allDesigns['StableStateAll'], allDesAng, onlysign = False)

stst = np.unique(allDesigns['StableStateAll'])
cmap2 = matl.cm.get_cmap('Set2',np.size(stst))
# cmap2 = matl.cm.get_cmap('jet',np.size(stst))
colors = cmap2(np.linspace(0,1,np.size(stst)))

#### Mask results out from the specified range of kappas
kappasrange = [10**-3,10**0]
dropers = ~((allDesigns.kappa < kappasrange[0]) | (allDesigns.kappa > kappasrange[-1]))
allDesigns = allDesigns[dropers]



#%%

symstst = ["" for x in stst]

for i in np.arange(np.size(stst,0)):
    allsym = allDesigns[allDesigns['StableStateAll'] == stst[i]]['Symetries'].unique()
    redsym = []
    for sym in allsym:
        redsym = np.append(redsym, sym[1:-1].split())
    symstst[i] = ' '.join(np.unique(redsym))
    print('Stable state',stst[i],'has symmetries:', symstst[i])
        
#%%
plt.close('all')   

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
        
    for j in stst:
        thisstst = thisang[thisang.StableStateAll == j]

        ax1.scatter(thisstst['kappa'], thisstst['face1'], c = colors[thisstst['StableStateAll']-1], s = 2)
        ax1.scatter(thisstst['kappa'], thisstst['face2'], c = colors[thisstst['StableStateAll']-1], s = 2)
        ax1.scatter(thisstst['kappa'], thisstst['face3'], c = colors[thisstst['StableStateAll']-1], s = 2)
        ax1.scatter(thisstst['kappa'], thisstst['face4'], c = colors[thisstst['StableStateAll']-1], s = 2)

    ax1.axhline(y=0.1, color='r', linestyle='-', linewidth = '0.4')
    
    fig1.show()
    fig1.savefig(Folder_name + '/Images/AreaFaces_' + str(restangles[i].astype(float))+'.pdf', transparent = True)
    fig1.savefig(Folder_name + '/Images/AreaFaces_' + str(restangles[i].astype(float))+ '.png', transparent = True)

#%%
plt.close('all')   

markers = np.array(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])

yticks = [[-1.2,0,1.0],[-38,-21,-9,0,10],[-2.4,-0.5,0]]
ylim = [[-1.5,1.5],[-40,13],[-2.7,0.2]]


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
    
#    ax1.set_yticks(np.linspace(0, maxTotEn,5))
#    ax1.set_ylim([-0.005, maxTotEn+0.005])
#    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2f'))
    ax1.set_yticks(yticks[i])
    ax1.set_ylim(ylim[i])
#    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.1f'))
    ax1.set_xscale('log')
    ax1.set_xticks([0.001,0.01,0.1,1])
    ax1.set_xlim([0.0007,1.5])
        
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
fig3 = plt.figure(figsize=(cm2inch(17.8), cm2inch(7)))
ax3 = plt.subplot(111,projection='3d')
ax3.set_xlim([-np.pi,np.pi])
ax3.set_ylim([-np.pi,np.pi])
ax3.set_zlim([-np.pi,np.pi])
    
order = [5,6,7,8]

markers = ['o','^','s']

for ang, i  in zip(restangles, markers):

    thisstate = allDesigns[allDesigns['restang'] == ang]
#    thisstate['desang0'] = np.zeros(np.shape(thisstate['desang1']))

    if not thisstate.empty:
        
        ax3.scatter(thisstate.iloc[:,order[0]].values,thisstate.iloc[:,order[1]].values,thisstate.iloc[:,order[2]].values, 
                    c = colors[thisstate['StableStateAll']-1], marker = i)
for i in stst:
    ax3.scatter(-400,-400,-400, c = [colors[i-1]], label = i)

plt.legend()

#%%
allDesigns[['kappa','Hinge Number','StableStateAll','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforAllImages.csv', index = False)
allDesigns.groupby('StableStateAll').apply(lambda df: df.sample(1, random_state = 0))[['kappa','Hinge Number','StableStateAll','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforStableStatesImages.csv', index = False)