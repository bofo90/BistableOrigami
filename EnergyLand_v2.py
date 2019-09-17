# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:04:51 2019

@author: iniguez
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matl
import xarray as xr
import configparser
import os.path

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
        
        
def countStableStates(finalAngles, distance, plot = False):
    
    import scipy.cluster.hierarchy as hierarch
    from scipy.spatial.distance import pdist
                                       
    Z = hierarch.linkage(finalAngles, 'centroid')
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

def ReadMetadata(file):
    
    if os.path.isfile(file):
        metadata = configparser.RawConfigParser()
        metadata.read(file)
    
        restang = float(metadata.get('options','restang'))
    else:
        raise FileNotFoundError('No metafile found at the given directory. Changes to the script to put manually the variables are needed\n') 

    return restang

kappas = np.logspace(-3,1,5)#23

Folder_name = "Results/SingleVertex3/sqp/energy/01-Aug-2019DesignAnalysis"
file_name1 = "/EnergyData.csv" 
file_name2 = "/Hinges.csv"
file_name3 = "/PosStad.csv"
file_name4 = "/Angles.csv"

for subdir in os.listdir(Folder_name):
    
    plt.close('all')
    prevFolder_name = Folder_name + '/' + subdir
    
    allData = pd.DataFrame()
#    allAngles = np.empty((0,4))
    allAngles = np.empty((0,3))
    allStSt = pd.DataFrame()
    
    for k in kappas:
        folder_name = prevFolder_name + "/kh%.5f_kta1000.00_ke1.00_kf100.00" %k
    
        dataEnergy = pd.read_csv(folder_name+file_name1)
        simLen = np.size(dataEnergy['Hinge Number'])
        
        dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
        dataEnergy['HingeEnergy'] = dataEnergy['HingeEnergy']/k
        dataEnergy['kappa'] = np.ones(simLen)*k
    
        dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
        dataAngles = np.delete(dataAngles, 0, 1)
        dataAnglesOrd = orderAngles(dataAngles, 2, simLen)
        dataEnergy['StableStates'] = np.zeros((simLen,1))
        dataEnergy['StableStates'] = countStableStates(dataAnglesOrd, 1)
#        dataEnergy[['ang1','ang2','ang3','ang4']] = pd.DataFrame(dataAnglesOrd)
        dataEnergy[['ang1','ang2','ang3']] = pd.DataFrame(dataAnglesOrd)
        
#        allStSt  = allStSt.append(dataEnergy.groupby('StableStates')[['ang1','ang2','ang3','ang4']].mean(), ignore_index = True)
        allStSt  = allStSt.append(dataEnergy.groupby('StableStates')[['ang1','ang2','ang3']].mean(), ignore_index = True)
        allAngles = np.append(allAngles, dataAngles, axis = 0)
        allData = allData.append(dataEnergy)
        
    restang = ReadMetadata(folder_name+'/metadata.txt')/np.pi   
     
    TotSimul = allData.shape[0]
    
    print(allData['Flags'].value_counts())
    exfl = allData.Flags.values.reshape(TotSimul,1)
    flagmask = (exfl !=1) & (exfl !=2)
    flagmask = ~flagmask.any(axis= 1)
    allData['Mask'] = flagmask
    
    allData.set_index(['kappa','Hinge Number'], inplace = True)
    allData.sort_index(inplace=True)
    
    kappasStSt = allData.groupby('kappa').apply(lambda _df: _df.groupby('StableStates')[['EdgeEnergy', 'HingeEnergy', 'TotalEnergy']].mean())
    kappasStSt = kappasStSt.reset_index(level=0)
    kappasStSt = kappasStSt.reset_index(level=0, drop=True)
#    kappasStSt['StableState'] = countStableStates(allStSt[['ang1','ang2','ang3','ang4']],1)
    kappasStSt['StableState'] = countStableStates(allStSt[['ang1','ang2','ang3']],1.1)
    
    
    selection = allData.groupby('kappa', as_index=False).apply(lambda _df: _df.groupby('StableStates').apply(lambda _df2: _df2.sample(1, random_state = 0)))
    selection = selection.reset_index(level = [0,1], drop = True)
    selection = selection.reset_index(level = [0,1])
    
    kappasStSt = kappasStSt.join(selection['Hinge Number'])
    kappasStSt['LogKappas'] = np.log10(kappasStSt.kappa)
    
    
    data = kappasStSt.to_xarray()
    
    
    
    #%%
    fig1 = plt.figure(1,figsize=(cm2inch(17.8), cm2inch(7)))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133)
    fig1.subplots_adjust(top=0.99,
    bottom=0.17,
    left=0.07,
    right=0.98,
    hspace=0.25,
    wspace=0.295)
    
    stst = np.unique(data.StableState)
    
    cmap = matl.cm.get_cmap('Set2',np.size(stst))
    
    for i, j in zip(stst, cmap(np.linspace(0,1,np.size(stst)))):
        onstst = np.array(data.StableState == i)
        ax1.scatter(data.kappa[onstst], data.TotalEnergy[onstst], c = [j], label = i)#
        ax2.scatter(data.kappa[onstst], data.HingeEnergy[onstst], c = [j])
        ax3.scatter(data.kappa[onstst], data.EdgeEnergy[onstst], c = [j])
    
    
    NiceGraph2D(ax1, 'Kappa', 'Total Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],0.6], divisions=[np.nan, 3], buffer=[0, 0.01])
    NiceGraph2D(ax2, 'Kappa', 'Normalized Hinge Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],40], buffer=[0, 0.05])
    NiceGraph2D(ax3, 'Kappa', 'Normalized Edge Energy', mincoord=[kappas[0],0], maxcoord=[kappas[-1],0.6], divisions=[np.nan, 3], buffer=[0, 0.01])
    
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax3.set_xscale('log')
    
    leg = ax1.legend(loc = 2, fontsize = 7, framealpha = 0.8, edgecolor = 'inherit', fancybox = False) 
    #           borderpad = 0.3, labelspacing = 0.1, handlelength = 0.4, handletextpad = 0.4)
    plt.setp(leg.get_texts(), color='0.2')
    leg.get_frame().set_linewidth(0.4)
    
    
    #%%
    fig1.show()
    fig1.savefig(prevFolder_name + '/Energy.pdf', transparent = True)
    fig1.savefig(prevFolder_name + '/Energy.png', transparent = True)
