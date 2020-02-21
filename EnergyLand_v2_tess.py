# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:04:51 2019

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
import ReadAndAnalyze as raa
matl.rcParams['pdf.fonttype'] = 42
matl.rcParams['ps.fonttype'] = 42
matl.rcParams['font.family'] = 'sans-serif'
matl.rcParams['font.sans-serif'] = 'Arial'
matl.rcParams['mathtext.fontset'] = 'cm'


def cm2inch(value):
    return value/2.54

def NiceGraph2D(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],buffer = [0.0, 0.0, 0.0]):
    
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

def NiceGraph2Dlog(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],buffer = [0.0, 0.0, 0.0]):
    
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
plt.close('all')

Folder_name = "Results/Tessellation4/25/2CFF/09-Jan-2020_5_5_"

allDesigns = pd.DataFrame()
# allKappasAnalysis = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    
tessellation = np.array(Folder_name.split('_')[1:-1]).astype(int)
numUC = np.prod(tessellation)
  
for subdir in os.listdir(Folder_name):
    if subdir == 'Images':
        continue    
        
    for subdir2 in os.listdir(Folder_name+'/'+subdir):
        folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
        
        ThisData, simLen = raa.ReadFile(folder_name)
        ThisData = raa.maskBadResults(ThisData)  
        
        ThisData.iloc[:,-numUC*4::] = raa.orderAnglesMat(ThisData.iloc[:,-numUC*4::].to_numpy(), 4, tessellation)
        ThisData['StableStates'] = raa.countStableStates(ThisData.iloc[:,-numUC*4::], 0.6, 'centroid')
        
        selData = raa.makeSelectionPerStSt(ThisData, simLen*0.05)
        
        allDesigns = allDesigns.append(selData)
        
       
    # kappasnumStSt = kappasStSt.groupby('kappa')['StableStates'].nunique()
    # kappasnumStSt = kappasnumStSt.reset_index()
    # kappasnumStSt['restang'] = np.ones(np.shape(kappasnumStSt)[0])*restang

    # allKappasAnalysis = allKappasAnalysis.append(kappasnumStSt)

allDesigns = allDesigns.reset_index(level=0, drop =True)
allDesigns = allDesigns.round(8)
restangles = allDesigns.restang.drop_duplicates().values
kappas = allDesigns.kappa.drop_duplicates().values


#%%
posStSt = np.shape(allDesigns)[0]
allUCang = np.resize(allDesigns.iloc[:,8:8+numUC*4].values,(posStSt*numUC,4))
[allUCang,sym] = raa.orderAngles(allUCang, 4, posStSt*numUC)
matUCStSt = raa.countStableStates(allUCang, 0.7, 'centroid',True)
matUCStSt = raa.standarizeStableStates(matUCStSt, allUCang, onlysign = False)
matUCStSt = np.reshape(matUCStSt, (posStSt, numUC))

allDesigns['StableStateAll'] = raa.countStableStates(allDesigns.iloc[:,5:5+numUC*4].values, 1, 'centroid')
[allDesigns['StableStateFromUC'],allDesigns['StableStatefUCName']] = raa.extractStableStates(matUCStSt)

stst = np.unique(allDesigns['StableStateAll'])
[ststfUC, nameloc] = np.unique(allDesigns['StableStateFromUC'], return_index = True)
ststfUCname = allDesigns['StableStatefUCName'][nameloc]
ststUC = np.unique(matUCStSt)

#### Mask results out from the specified range of kappas
# kappasrange = [10**-3,10**0]
# dropers = ~((allDesigns.kappa < kappasrange[0]) | (allDesigns.kappa > kappasrange[-1]))
# allDesigns = allDesigns[dropers]

# cmap2 = matl.cm.get_cmap('Set2',np.size(stst))
cmap2 = matl.cm.get_cmap('jet',np.size(stst))
colors = cmap2(np.linspace(0,1,np.size(stst)))

cmapfUC = matl.cm.get_cmap('Set3',np.size(ststfUC))
# cmapfUC = matl.cm.get_cmap('jet',np.size(ststfUC))
colorsfUC = cmapfUC(np.linspace(0,1,np.size(ststfUC)))

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
    NiceGraph2D(ax1, r'$\kappa$', r'$K$')
    
# #    ax1.set_yticks(np.linspace(0, maxTotEn,5))
# #    ax1.set_ylim([-0.005, maxTotEn+0.005])
# #    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2f'))
#     ax1.set_yticks(yticks[i])
#     ax1.set_ylim(ylim[i])
# #    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.1f'))
    ax1.set_xscale('log')
    ax1.set_xticks([0.001,0.01,0.1,1])
    ax1.set_xlim([0.0007,1.5])
        
#     for j in stst:
#         thisstst = thisang[thisang.StableStateAll == j]
# #        ax1.scatter(thisstst['kappa'], thisstst['TotalEnergy'], c = colors[thisstst['StableStateAll'].values.astype('int')-1], s = 10)
#         ax1.scatter(thisstst['kappa'], thisstst['Curvature'], c = colors[thisstst['StableStateAll']-1], s = 10)
        
    for j in ststfUC:
        thisstst = thisang[thisang.StableStateFromUC == j]
#        ax1.scatter(thisstst['kappa'], thisstst['TotalEnergy'], c = colorsfUC[thisstst['StableStateFromUC'].values.astype('int')-1], s = 10)
        ax1.scatter(thisstst['kappa'], thisstst['Curvature'], c = colorsfUC[thisstst['StableStateFromUC']-1], s = 10)
 
    fig1.show()
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+'.pdf', transparent = True)
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+ '.png', transparent = True)
    fig1.savefig(Folder_name + '/Images/Restang_' + str(restangles[i].astype(float))+'.pdf', transparent = True)
    fig1.savefig(Folder_name + '/Images/Restang_' + str(restangles[i].astype(float))+ '.png', transparent = True)


cmapfig, normfig = from_levels_and_colors(np.arange(np.size(ststfUC)+1),matl.cm.Set3(np.linspace(0, 1,np.size(ststfUC))))

fig1b = plt.figure(figsize=(cm2inch(8.8), cm2inch(1.5)))
fig1b.subplots_adjust(top=0.884,
bottom=0.116,
left=0.039,
right=0.961)
ax1b = plt.subplot(111)
cbar2 = plt.colorbar(matl.cm.ScalarMappable(norm=normfig, cmap=cmapfig), ax = ax1b, fraction=0.99, pad=0.01, orientation='horizontal')#, extend = 'max'
cbar2.set_ticks(ststfUC-0.5)
cbar2.ax.set_xticklabels(ststfUCname.astype(str))
cbar2.set_label('Stable State', fontsize = 9, color = '0.2',labelpad = 0)
cbar2.ax.tick_params(colors='0.2', pad=2)
cbar2.outline.set_edgecolor('0.2')
cbar2.outline.set_linewidth(0.4)
ax1b.remove()   

fig1b.savefig(Folder_name + '/Images/Restang_CB.pdf', transparent = True)
fig1b.savefig(Folder_name + '/Images/Restang_CB.png', transparent = True)


#%%
fig3 = plt.figure(figsize=(cm2inch(10), cm2inch(7)))
ax3 = plt.subplot(111,projection='3d')
ax3.set_xlim([-np.pi,np.pi])
ax3.set_ylim([-np.pi,np.pi])
ax3.set_zlim([-np.pi,np.pi])

markers = np.array(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
    
ax3.scatter(allUCang[:,0],allUCang[:,1],allUCang[:,2], c = colorsUC[matUCStSt.flatten()-1])
for i in ststUC:
    ax3.scatter(-400,-400,-400, c = [colorsUC[i-1]], label = i)

# for i, j in zip(np.arange(numUC), markers):
#     ax3.scatter(allUCang[i::numUC,0],allUCang[i::numUC,1],allUCang[i::numUC,2], 
#                 c = colors[allDesigns['StableStateAll'].values-1], marker = j)
# for i in stst:
#     ax3.scatter(-400,-400,-400, c = [colors[i-1]], label = i)
    
# for i in np.arange(numUC):
#     ax3.scatter(allUCang[i::numUC,0],allUCang[i::numUC,1],allUCang[i::numUC,2], 
#                 c = colorsfUC[allDesigns['StableStateFromUC'].values-1])
# for i in ststfUC:
#     ax3.scatter(-400,-400,-400, c = [colorsfUC[i-1]], label = i)

plt.legend()

#%%
allDesigns[['kappa','Hinge Number','StableStatefUCName','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforAllImages.csv', index = False)
allDesigns.groupby('StableStateAll').apply(lambda df: df.sample(1, random_state = 0))[['kappa','Hinge Number','StableStatefUCName','restang','Curvature']].to_csv(Folder_name + '/Images/InfoforStableStatesImages.csv', index = False)