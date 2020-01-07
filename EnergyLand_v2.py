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
import xarray as xr
import configparser
import os.path
import ternary #from https://github.com/marcharper/python-ternary
import itertools
import copy
matl.rcParams['pdf.fonttype'] = 42
matl.rcParams['ps.fonttype'] = 42
matl.rcParams['font.family'] = 'sans-serif'
matl.rcParams['font.sans-serif'] = 'Arial'
matl.rcParams['mathtext.fontset'] = 'cm'


def cm2inch(value):
    return value/2.54

def orderAngles(angles, steps, simulations):
    
    symetries = ["" for x in range(simulations)]
    finalAngles = np.zeros((simulations,np.size(angles,1)))
    for hinge in np.arange(simulations):
        sym = ''
        if np.sum(np.sign(np.around(angles[hinge,:], decimals = 4))) < 0: ### sort to have mayority positive angles
            angles[hinge,:] = angles[hinge,::-1]*(-1)
            sym = sym + 'm'
        angles[hinge,:]= np.roll(angles[hinge,:],np.argmin(angles[hinge,:]))
        if np.argmin(angles[hinge,:]) != 0:
            sym = sym+'r'+ str(np.argmin(angles[hinge,:]))
        # angles[hinge,:] = np.sort(angles[hinge,:])          ### sort from smallest to highest angles
        finalAngles[hinge,:] = angles[hinge,:]
        symetries[hinge] = sym
        
    return [finalAngles, symetries]
        
        
def countStableStates(finalAngles, distance, method, plot = False):
    
    import scipy.cluster.hierarchy as hierarch
    from scipy.spatial.distance import pdist
                                       
    Z = hierarch.linkage(finalAngles, method)
    inverse = hierarch.fcluster(Z, distance, criterion='distance')
    c = hierarch.cophenet(Z, pdist(finalAngles))

    if plot:
        
        print('this is the cophenet of the hierarchical linkage', c[0])
        
        plt.figure(0,figsize=(25, 10))
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
    axes.tick_params(axis='x',colors=gray, direction = 'in', width = 0.4)
    axes.spines['bottom'].set_color(gray)
    axes.spines['top'].set_color(gray)
    axes.yaxis.label.set_color(gray)
    axes.tick_params(axis='y', colors=gray, direction = 'in', width = 0.4)
    axes.spines['left'].set_color(gray)
    axes.spines['right'].set_color(gray)
    axes.tick_params(pad = 2)
    
    axes.tick_params(axis='y', which='minor', colors='0.2', direction = 'in', width = 0.4)
    axes.tick_params(axis='x', which='minor', colors='0.2', direction = 'in', width = 0.4)

    
    for axis in ['top','bottom','left','right']:
        axes.spines[axis].set_linewidth(0.4)
        
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
def NiceTerciaryGraph(ax, name, scale, divisions):
    
    ax.axis('off')
    tax = ternary.TernaryAxesSubplot(ax=ax, scale = scale)
    
    fontsize = 9
    matl.rcParams.update({'font.size': 9})
    gray = '0.2'
    linewidth = 0.4
    mult = scale/divisions
    
    tax.boundary(linewidth=linewidth)
    tax.gridlines(color=gray, multiple=mult)
    
    tax.set_title(name, fontsize=fontsize, color = gray, pad = 15)
    tax.left_axis_label("Angle1", fontsize=fontsize, color = gray, offset = 0.3)
    tax.right_axis_label("Angle2", fontsize=fontsize, color = gray, offset = 0.3)
    tax.bottom_axis_label("Angle3", fontsize=fontsize, color = gray, offset = 0.3)
    
    tax.ticks(axis='lbr', linewidth=linewidth, axes_colors = {'l': gray, 'r':gray, 'b': gray},fontsize = fontsize,
              ticks = [180,150,120,90,60,30,0], clockwise = True, offset = 0.05)
    tax.clear_matplotlib_ticks()
    
    return tax

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

Folder_name = "Results/SingleVertex4/2CFF/02-Dec-2019_0.00_90.00_180.00_270.00_"
file_name1 = "/EnergyData.csv" 
file_name2 = "/Hinges.csv"
file_name3 = "/PosStad.csv"
file_name4 = "/Angles.csv"

allDesigns = pd.DataFrame()
allKappasAnalysis = pd.DataFrame()

if not os.path.isdir(Folder_name + '/Images/'):
    os.mkdir(Folder_name + '/Images/')
    os.mkdir(Folder_name + '/Images/SingleDes/')
    
designang = np.array(Folder_name.split('_')[1:-1]).astype(float)
  
for subdir in os.listdir(Folder_name):
    if subdir == 'Images':
        continue    
    
    allData = pd.DataFrame()
    allAngles = np.empty((0,4))
#    allAngles = np.empty((0,3))    
    
    for subdir2 in os.listdir(Folder_name+'/'+subdir):

        folder_name = Folder_name+'/'+subdir+'/'+subdir2+'/energy'
    
        dataEnergy = pd.read_csv(folder_name+file_name1)
        simLen = np.size(dataEnergy['Hinge Number'])
        
        dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
        dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']
        dataEnergy['kappa'] = np.ones(simLen)*np.float(subdir2.split('_')[1])
        
        dataHinges = pd.read_csv(folder_name+file_name2)
        dataEnergy['Curvature'] = dataHinges['Curvature']
    
        dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
        dataAngles = np.delete(dataAngles, 0, 1)
        [dataAnglesOrd, anglesSym] = orderAngles(dataAngles, 2, simLen)
        dataEnergy['StableStates'] = np.zeros((simLen,1))
        dataEnergy['StableStates'] = countStableStates(dataAnglesOrd, 10, 'ward')
#        dataEnergy['StableStates'] = countStableStates(dataAnglesOrd, 0.5, 'centroid')
        dataEnergy[['ang1','ang2','ang3','ang4']] = pd.DataFrame(dataAnglesOrd)
#        dataEnergy[['ang1','ang2','ang3']] = pd.DataFrame(dataAnglesOrd)
        dataEnergy['Symmetry'] = anglesSym
        
        allAngles = np.append(allAngles, dataAngles, axis = 0)
        allData = allData.append(dataEnergy)
        
    restang = np.float(subdir.split('_')[1])
    
     
    TotSimul = allData.shape[0]
    
    print(allData['Flags'].value_counts())
    exfl = allData.Flags.values.reshape(TotSimul,1)
    flagmask = (exfl !=1) & (exfl !=2)
    flagmask = ~flagmask.any(axis= 1)
    allData['Mask'] = flagmask
    
    allData.set_index(['kappa','Hinge Number'], inplace = True)
    allData.sort_index(inplace=True)
    
    kappasStSt = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['EdgeEnergy', 'HingeEnergy', 'TotalEnergy','ang1','ang2','ang3','ang4', 'Curvature']].mean())
#    kappasStSt = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['EdgeEnergy', 'HingeEnergy', 'TotalEnergy','ang1','ang2','ang3', 'Curvature']].mean())
    kappasStSt['amountSym'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['Symmetry']].nunique())
    kappasStSt['Symetries'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates').apply(lambda x: str(np.sort(x.Symmetry.unique()))))
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['DiagonalEnergy']].count())
    more10percent = kappasStSt['amountStSt']>simLen*0.1
    kappasStSt = kappasStSt[more10percent]
    
    kappasStSt = kappasStSt.reset_index(level=0)
    kappasStSt = kappasStSt.reset_index(level=0)#, drop=True)
#    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3','ang4']],1, 'centroid')
#    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3']],0.45, 'centroid')
    
    selection = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates').apply(lambda _df2: _df2.sample(1, random_state = 0)))
    selection = selection[more10percent]
    selection = selection.reset_index(level = [0,1], drop = True)
    selection = selection.reset_index(level = [0,1])
    
    kappasStSt = kappasStSt.join(selection['Hinge Number'])
    kappasStSt['LogKappas'] = np.log10(kappasStSt.kappa)
    
    kappasStSt['restang'] = np.ones(np.shape(kappasStSt)[0])*restang
    
    kappasnumStSt = kappasStSt.groupby('kappa')['StableStates'].nunique()
    kappasnumStSt = kappasnumStSt.reset_index()
    kappasnumStSt['restang'] = np.ones(np.shape(kappasnumStSt)[0])*restang

    allKappasAnalysis = allKappasAnalysis.append(kappasnumStSt)
    
    data = kappasStSt.to_xarray()
    
    allDesigns = allDesigns.append(kappasStSt)    
      
allDesigns = allDesigns.reset_index(level=0, drop =True)
restangles = allKappasAnalysis.restang.drop_duplicates().values


#%%
allDesigns['StableStateAll'] = countStableStates(allDesigns[['ang1','ang2','ang3','ang4']], 0.8, 'centroid', True)
#allDesigns['StableStateAll'] = countStableStates(allDesigns[['ang1','ang2','ang3']], 0.5, 'centroid')
stst = np.unique(allDesigns['StableStateAll'])
############################################# TO make SS consistence between plots (its done manualy)
# newSS = np.array([4,4,3,3,1,1,2,1,2])
# invmask_copy = copy.deepcopy(allDesigns.StableStateAll)

# for old,new in zip(stst, newSS):
#     allDesigns.StableStateAll[invmask_copy == old] = new

# delSS = np.array([4,5,6,7,8,9])
# stst = np.delete(stst, delSS)
#####################################################################################################


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

markers = np.array(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
cmap = matl.cm.get_cmap('jet',np.size(stst))
colors = cmap(np.linspace(0,1,np.size(stst)))

yticks = [[-1.2,-0.6,0,0.6,1.2],[-36,-24,-12,0,12],[-2.4,-1.6,-0.8,0]]
ylim = [[-1.5,1.5],[-39,13],[-2.7,0.2]]


for i in np.arange(np.size(restangles)):
    
    fig1 = plt.figure(figsize=(cm2inch(5.5), cm2inch(4)))
    ax1 = plt.subplot(111)
    fig1.subplots_adjust(top=0.962,
    bottom=0.25,
    left=0.22,
    right=0.987)
    
    thisangBool = allDesigns['restang'] == restangles[i]
    thisang = allDesigns[thisangBool]    
    
    maxTotEn = thisang['TotalEnergy'].max()
    maxCurv = np.ceil(thisang['Curvature'].max())#quantile(0.95)
    minCurv = np.floor(thisang['Curvature'].min())
    
#    NiceGraph2D(ax1, r'$\kappa$', r'$E_\mathrm{tot}$')
    NiceGraph2D(ax1, r'$\kappa$', r'$C$')#, mincoord=[np.log10(kappas[0]),-3], maxcoord=[np.log10(kappas[-1]),3], divisions=[5, 5], buffer=[0.1, 1])
    
#    ax1.set_yticks(np.linspace(0, maxTotEn,5))
#    ax1.set_ylim([-0.005, maxTotEn+0.005])
#    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2f'))
    ax1.set_yticks(yticks[i])
    ax1.set_ylim(ylim[i])
#    ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.1f'))
    ax1.set_xscale('log')
    ax1.set_xticks([0.001,0.01,0.1,1,10])
    ax1.set_xlim([0.0007,15])
        
    for j in stst:
        thisstst = thisang[thisang.StableStateAll == j]
#    plt.close('all')

#        ax1.scatter(thisstst['kappa'], thisstst['TotalEnergy'], c = colors[thisstst['StableStateAll'].values.astype('int')-1], s = 10)
        ax1.scatter(thisstst['kappa'], thisstst['Curvature'], c = colors[thisstst['StableStateAll'].values.astype('int')-1], s = 5)#,
#                s = 5, marker = markers[i])#, linestyle = lines[i], lw = 2.5)

#leg = ax1.legend(loc = 2, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False) 
##           borderpad = 0.3, labelspacing = 0.1, handlelength = 0.4, handletextpad = 0.4)
#plt.setp(leg.get_texts(), color='0.2')
#leg.get_frame().set_linewidth(0.4)
    
    fig1.show()
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+'.pdf', transparent = True)
#    fig1.savefig(Folder_name + '/Images/Energy_Restang_' + str(i.astype(float))+ '.png', transparent = True)
    # fig1.savefig(Folder_name + '/Images/Restang_' + str(restangles[i].astype(float))+'.pdf', transparent = True)
    # fig1.savefig(Folder_name + '/Images/Restang_' + str(restangles[i].astype(float))+ '.png', transparent = True)

#%%
fig3 = plt.figure(figsize=(cm2inch(17.8), cm2inch(7)))
ax3 = plt.subplot(111,projection='3d')
ax3.set_xlim([-np.pi,np.pi])
ax3.set_ylim([-np.pi,np.pi])
ax3.set_zlim([-np.pi,np.pi])
    
#cmap2 = matl.cm.get_cmap('Set2',np.size(kappas))
#colors = cmap2(np.linspace(0,1,np.size(kappas)))    
cmap2 = matl.cm.get_cmap('jet',np.size(stst))
colors = cmap2(np.linspace(0,1,np.size(stst)))

order = [5,6,7,8]

markers = ['o','^','s']

for ang, i  in zip(restangles, markers):

    thisstate = allDesigns[allDesigns['restang'] == ang]
#    thisstate['desang0'] = np.zeros(np.shape(thisstate['desang1']))

    if not thisstate.empty:
        
        ax3.scatter(thisstate.iloc[:,order[0]].values,thisstate.iloc[:,order[1]].values,thisstate.iloc[:,order[2]].values, 
                    c = colors[thisstate['StableStateAll']-1], marker = i)

#%%
allDesigns[['kappa','Hinge Number','StableStateAll','restang']].to_csv(Folder_name + '/Images/InfoforAllImages.csv', index = False)
allDesigns.groupby('StableStateAll').apply(lambda df: df.sample(1))[['kappa','Hinge Number','StableStateAll','restang']].to_csv(Folder_name + '/Images/InfoforStableStatesImages.csv', index = False)