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
import ternary

def cm2inch(value):
    return value/2.54

def orderAngles(angles, steps, simulations):
    
    finalAngles = np.empty((0,np.size(angles,1)))
    for hinge in np.arange(simulations):
        if np.sum(np.sign(np.around(angles[steps*(hinge+1)-1,:], decimals = 4))) < 0:
            angles[steps*(hinge+1)-1,:] = angles[steps*(hinge+1)-1,:]*(-1)        
        sortAllAngIndex = np.lexsort((angles[steps*(hinge+1)-1,:],angles[steps*hinge,:]))#[0,1,2,3]
        finalAngles = np.append(finalAngles, [angles[steps*(hinge+1)-1,sortAllAngIndex]], axis = 0)
        
    return finalAngles
        
        
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

def NiceGraph2Dlog(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
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
        axes.set_ylim([mincoord[1]*((buffer[1])**(-1)), maxcoord[1]*(buffer[1])])
        if isinstance(divisions[1], (list, tuple, np.ndarray)):
            if ~np.isnan(divisions[1]).any():
                axes.set_yticks(divisions[1])
        else:
            if ~np.isnan(divisions[1]):
                axes.set_yticks(np.logspace(np.log10(mincoord[1]),np.log10(maxcoord[1]),divisions[1]))
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
        metadata = configparser.RawConfigParser()
        metadata.read(file)
    
        restang = float(metadata.get('options','restang'))
        designang = np.array(metadata.get('options','angDesign').split(),dtype = int)
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

    return restang/np.pi, designang

kappas = np.logspace(-3,1,5)#23

Folder_name = "Results/SingleVertex3/sqp/energy/01-Aug-2019DesignAnalysis"
file_name1 = "/EnergyData.csv" 
file_name2 = "/Hinges.csv"
file_name3 = "/PosStad.csv"
file_name4 = "/Angles.csv"

allDesigns = pd.DataFrame()

for subdir in os.listdir(Folder_name):
    
    if subdir == 'Images':
        continue
    
    plt.close('all')
    prevFolder_name = Folder_name + '/' + subdir
    
    allData = pd.DataFrame()
#    allAngles = np.empty((0,4))
    allAngles = np.empty((0,3))
    
    for k in kappas:
        folder_name = prevFolder_name + "/kh%.5f_kta1000.00_ke1.00_kf100.00" %k
    
        dataEnergy = pd.read_csv(folder_name+file_name1)
        simLen = np.size(dataEnergy['Hinge Number'])
        
        dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
        dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']
        dataEnergy['kappa'] = np.ones(simLen)*k
    
        dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
        dataAngles = np.delete(dataAngles, 0, 1)
        dataAnglesOrd = orderAngles(dataAngles, 2, simLen)
        dataEnergy['StableStates'] = np.zeros((simLen,1))
        dataEnergy['StableStates'] = countStableStates(dataAnglesOrd, 0.3, 'centroid')
#        dataEnergy[['ang1','ang2','ang3','ang4']] = pd.DataFrame(dataAnglesOrd)
        dataEnergy[['ang1','ang2','ang3']] = pd.DataFrame(dataAnglesOrd)
        
        allAngles = np.append(allAngles, dataAngles, axis = 0)
        allData = allData.append(dataEnergy)
        
    restang, designang = ReadMetadata(folder_name+'/metadata.txt') 
     
    TotSimul = allData.shape[0]
    
    print(allData['Flags'].value_counts())
    exfl = allData.Flags.values.reshape(TotSimul,1)
    flagmask = (exfl !=1) & (exfl !=2)
    flagmask = ~flagmask.any(axis= 1)
    allData['Mask'] = flagmask
    
    allData.set_index(['kappa','Hinge Number'], inplace = True)
    allData.sort_index(inplace=True)
    
    kappasStSt = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['EdgeEnergy', 'HingeEnergy', 'TotalEnergy','ang1','ang2','ang3']].mean())
    ####Only select stable states that have more than 10% appearance in each simulation
    kappasStSt['amountStSt'] = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['DiagonalEnergy']].count())
    kappasStSt = kappasStSt[kappasStSt['amountStSt']>simLen*0.1]
    
    kappasStSt = kappasStSt.reset_index(level=0)
    kappasStSt = kappasStSt.reset_index(level=0, drop=True)
#    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3','ang4']],1, 'centroid')
    kappasStSt['StableState'] = countStableStates(kappasStSt[['ang1','ang2','ang3']],0.45, 'centroid')
    
    selection = allData.groupby('kappa', as_index=False).apply(lambda _df: _df.groupby('StableStates').apply(lambda _df2: _df2.sample(1, random_state = 0)))
    selection = selection.reset_index(level = [0,1], drop = True)
    selection = selection.reset_index(level = [0,1])
    
    kappasStSt = kappasStSt.join(selection['Hinge Number'])
    kappasStSt['LogKappas'] = np.log10(kappasStSt.kappa)
    
    kappasStSt['desang1'] = np.ones(np.shape(kappasStSt)[0])*designang[1]
    kappasStSt['desang2'] = np.ones(np.shape(kappasStSt)[0])*(designang[2]-designang[1])
    
    data = kappasStSt.to_xarray()
    
    allDesigns = allDesigns.append(kappasStSt)    
    
    #%%
#    fig1 = plt.figure(1,figsize=(cm2inch(17.8), cm2inch(7)))
#    ax1 = plt.subplot(131)
#    ax2 = plt.subplot(132)
#    ax3 = plt.subplot(133)
#    fig1.subplots_adjust(top=0.99,
#    bottom=0.17,
#    left=0.07,
#    right=0.98,
#    hspace=0.25,
#    wspace=0.295)
#    
#    stst = np.unique(data.StableState)
#    
#    cmap = matl.cm.get_cmap('Set2',np.size(stst))
#    
#    for i, j in zip(stst, cmap(np.linspace(0,1,np.size(stst)))):
#        onstst = np.array(data.StableState == i)
#        ax1.scatter(data.kappa[onstst], data.TotalEnergy[onstst], c = [j], label = i)#
#        ax2.scatter(data.kappa[onstst], data.HingeEnergy[onstst], c = [j])
#        ax3.scatter(data.kappa[onstst], data.EdgeEnergy[onstst], c = [j])
#    
#    
#    NiceGraph2D(ax1, 'Kappa', 'Total Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],1], divisions=[np.nan, 3], buffer=[0, 0.01])
#    NiceGraph2D(ax2, 'Kappa', 'Normalized Hinge Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],1], buffer=[0, 0.05])
#    NiceGraph2D(ax3, 'Kappa', 'Normalized Edge Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],1], divisions=[np.nan, 3], buffer=[0, 0.01])
#    
#    ax1.set_xscale('log')
#    ax2.set_xscale('log')
#    ax3.set_xscale('log')
#    
#    leg = ax1.legend(loc = 2, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False) 
#    #           borderpad = 0.3, labelspacing = 0.1, handlelength = 0.4, handletextpad = 0.4)
#    plt.setp(leg.get_texts(), color='0.2')
#    leg.get_frame().set_linewidth(0.4)
#    
#    
#    #%%
#    fig1.show()
#    fig1.savefig(Folder_name + '/Images/' + subdir[7:] + '.pdf', transparent = True)
#    fig1.savefig(Folder_name + '/Images/' + subdir[7:] + '.png', transparent = True)
    
#%%
allDesigns['desang3'] = 360-allDesigns['desang1']-allDesigns['desang2']
allDesigns = allDesigns.reset_index(level=0, drop =True)

#%%
allDesigns['StableStateAll'] = countStableStates(allDesigns[['ang1','ang2','ang3']], 4, 'ward')
stst = np.unique(allDesigns['StableStateAll'])
cmap2 = matl.cm.get_cmap('Set2',np.size(stst))
colors = cmap2(np.linspace(0,1,np.size(stst)))

for state in stst:
    plt.close('all')
    
    fig2, axes = plt.subplots(2,3, figsize=(cm2inch(15.8), cm2inch(10)))
    fig2.subplots_adjust(top=0.91,
    bottom=0.055,
    left=0.035,
    right=0.97,
    hspace=0.41,
    wspace=0.26)

    thisstate = allDesigns[allDesigns['StableStateAll'] == state]
    
    for i, ax in enumerate(axes.flat):
        if i >= np.size(kappas):
            ax.axis('off')
            continue
        thiskappa = thisstate[thisstate['kappa'] == kappas[i]]
        
        tax = NiceTerciaryGraph(ax, 'kappa '+str(kappas[i]) , 180, 6)
        
        if not thiskappa.empty:
            tax.scatter(-thiskappa[['desang1','desang2','desang3']].values+180, c = colors[thiskappa['StableStateAll']-1], s = 4)
            tax.scatter(-thiskappa[['desang2','desang3','desang1']].values+180, c = colors[thiskappa['StableStateAll']-1], s = 4)
            tax.scatter(-thiskappa[['desang3','desang1','desang2']].values+180, c = colors[thiskappa['StableStateAll']-1], s = 4)
    
    fig2.show()
    fig2.savefig(Folder_name + '/Images/' + 'DesignSpaceStSt' + str(state) + '.pdf', transparent = True)
    fig2.savefig(Folder_name + '/Images/' + 'DesignSpaceStSt' + str(state) + '.png', transparent = True)

fig3 = plt.figure(3,figsize=(cm2inch(17.8), cm2inch(7)))
ax3 = plt.subplot(111,projection='3d')
for state in stst:
    thisstate = allDesigns[allDesigns['StableStateAll'] == state]
    ax3.scatter(thisstate['ang1'].values/np.pi,thisstate['ang2'].values/np.pi,thisstate['ang3'].values/np.pi, c = colors[thisstate['StableStateAll']-1])

allDesigns[['kappa','Hinge Number','desang1','desang2','desang3', 'StableStateAll']].to_csv(Folder_name + '/Images/InfoforImages.csv', index = False)
