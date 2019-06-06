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
        sortAllAngIndex = np.lexsort((angles[steps*(hinge+1)-1,:],angles[steps*hinge,:]))
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

kappas = np.unique(np.concatenate((np.logspace(-4,0,9), np.logspace(-4,0,9)*10**0.25)))


for k in kappas:
    plt.close('all')
    
    folder_name = "Results/SingleVertex4/sqp/energy/05-Jun-2019_kappas_analysis/kh%.5f_kta100.00_ke1.00_kf100.00" %k
    
    file_name1 = "/EnergyData.csv" 
    file_name2 = "/Hinges.csv"
    file_name3 = "/PosStad.csv"
    file_name4 = "/Angles.csv"
    
    dataEnergy = pd.read_csv(folder_name+file_name1)
    dataEnergy['TotalEnergy'] = dataEnergy['EdgeEnergy']+dataEnergy['DiagonalEnergy']+dataEnergy['HingeEnergy']
    
    dataVar = pd.read_csv(folder_name+file_name2)
    
    dataAngles = np.loadtxt(folder_name+file_name4,skiprows=1, delimiter = ',', dtype = np.float64)
    dataAngles = np.delete(dataAngles, 0, 1)
    
    TotSimul = dataVar.shape[0]
    IterPerSimulEnergy = np.int(dataEnergy.shape[0]/TotSimul)
    IterPerSimulAngles = np.int(dataAngles.shape[0]/TotSimul)
    
    print(dataEnergy['Flags'].value_counts())
    exfl = dataEnergy.Flags.values.reshape(TotSimul,IterPerSimulEnergy)
    flagmask = (exfl !=1) & (exfl !=2)
    flagmask = ~flagmask.any(axis= 1)
    
    dataAnglesOrd = orderAngles(dataAngles, IterPerSimulAngles, TotSimul)
    
    dataVar['Mask'] = flagmask
    dataVar['TargetAngle1'] = dataVar['TargetAngle1']/np.pi
    dataVar['TargetAngle2'] = dataVar['TargetAngle2']/np.pi
    dataVar['StableStates'] = np.zeros(np.size(dataVar['Mask']))
    dataVar['StableStates'] =  countStableStates(dataAnglesOrd, 1, True)
    
    np.random.seed(0)
    print(dataVar.groupby('StableStates')[['NumAngle1','NumAngle2']].apply(lambda x: x.sample(1,random_state = 0)))
    
    dataVar = dataVar.join(dataEnergy[['Hinge Number','TotalEnergy']][IterPerSimulEnergy-2::IterPerSimulEnergy].set_index('Hinge Number'), on = 'HingeNumber')
    dataVar = dataVar.join(dataEnergy[['Hinge Number','TargetAngleEnergy']][IterPerSimulEnergy-2::IterPerSimulEnergy].set_index('Hinge Number'), on = 'HingeNumber')
    dataVar.set_index(['TargetAngle1','TargetAngle2'], inplace = True)
    dataVar.sort_index(inplace=True)
    
    data = dataVar.to_xarray()
    
    #for i in data.Kappa.data:
    #    data['StableStates'][data.Kappa == i] = countStableStates(dataAnglesOrd[data['HingeNumber'][data.Kappa == i].data[0]-1,:], plot = False)
    #dataVar['StableStates'] = data.StableStates.data.flatten()
    
    
    #%%
    fig1 = plt.figure(1,figsize=(cm2inch(17.8), cm2inch(15)))
    ax1 = plt.subplot(111)
    fig1.subplots_adjust(top=0.995,
    bottom=0.07,
    left=0.1,
    right=0.915)
    
    minEnTA = np.float(data.TargetAngleEnergy.where(data.Mask).min())
    aveEnTA = np.float(data.TargetAngleEnergy.where(data.Mask).mean())
    stdEnTA = np.float(data.TargetAngleEnergy.where(data.Mask).std())
    
    energiesTA = np.linspace(minEnTA, aveEnTA + stdEnTA*2, 100)
    
    norm = matl.colors.LogNorm(10**energiesTA[0], 10**energiesTA[-1])
    cmap = matl.cm.ScalarMappable(norm = norm, cmap=matl.cm.summer)
    cmap.set_array([])
    
    data['TargetAngleEnergy'].where(data.Mask).plot(axes = ax1, cmap = matl.cm.summer, add_colorbar=True,
        levels = energiesTA)
    
    plt.title('')
    
    NiceGraph2D(ax1, 'Target Angle1', 'Target Angle2', mincoord = [-1,-1], maxcoord = [1,1],
                divisions = [5, 5], buffer = [0.01, 0.01])
    
    #cbar = fig1.colorbar(cmap, ticks = np.logspace(-4,0,5), fraction=0.05, pad=0.01) 
    #cbar.set_label('Kappa', fontsize = 9, color = '0.2')
    #cbar.ax.tick_params(axis='y',colors='0.2')
    #cbar.ax.tick_params(axis='x',colors='0.2')
    #cbar.outline.set_edgecolor('0.2')
    
    #%%
    fig2 = plt.figure(2,figsize=(cm2inch(17.8), cm2inch(15)))
    ax2 = plt.subplot(111)
    fig2.subplots_adjust(top=0.995,
    bottom=0.07,
    left=0.1,
    right=0.915)
    
    minEn = np.float(data.TotalEnergy.where(data.Mask).where(data.TargetAngleEnergy<0.5).min())
    aveEn = np.float(data.TotalEnergy.where(data.Mask).where(data.TargetAngleEnergy<0.5).mean())
    stdEn = np.float(data.TotalEnergy.where(data.Mask).where(data.TargetAngleEnergy<0.5).std())
    
    energies = np.log10(np.linspace(minEn,aveEn+stdEn*2,100))
    
    norm2 = matl.colors.LogNorm(10**energies[0], 10**energies[-1])
    cmap2 = matl.cm.ScalarMappable(norm = norm2, cmap=matl.cm.nipy_spectral)
    cmap2.set_array([])
    
    np.log10(data['TotalEnergy'].where(data.Mask)).where(data.TargetAngleEnergy<0.5).plot(axes = ax2, 
            cmap=matl.cm.nipy_spectral ,add_colorbar=False, levels = energies)
    
    NiceGraph2D(ax2, 'Target Angle1', 'Target Angle2', mincoord = [-1,-1], maxcoord = [1,1],
                divisions = [5, 5], buffer = [0.01, 0.01])
    
    cbar2 = fig2.colorbar(cmap2,  fraction=0.05, pad=0.01, extend = 'both')# ticks = np.logspace(-1,2,10),
    cbar2.set_label('Energy', fontsize = 9, color = '0.2')
    cbar2.ax.tick_params(axis='y',colors='0.2')
    cbar2.ax.tick_params(axis='x',colors='0.2')
    cbar2.outline.set_edgecolor('0.2')
    
    #%%
    fig3 = plt.figure(3,figsize=(cm2inch(17.8), cm2inch(15)))
    ax3 = plt.subplot(111)
    fig3.subplots_adjust(top=0.98,
    bottom=0.07,
    left=0.085,
    right=0.91)
    
    data['StableStates'].where(data.Mask).where(data.TargetAngleEnergy<0.5).plot(axes = ax3, cmap=matl.cm.Set2, add_colorbar=True,
        levels = np.linspace(data.StableStates.min()-0.5,data.StableStates.max()+0.5,int(data.StableStates.max())+1), robust = True,
        cbar_kwargs={'ticks': np.linspace(data.StableStates.min(),data.StableStates.max(),int(data.StableStates.max()))})
    
    NiceGraph2D(ax3, 'Target Angle1', 'Target Angle2', mincoord = [-1,-1], maxcoord = [1,1],
                divisions = [5, 5], buffer = [0.01, 0.01])
    
    #%%
    fig1.show()
    fig1.savefig(folder_name + '/Energy_Kappas.pdf', transparent = True)
    fig1.savefig(folder_name + '/Energy_Kappas.png', transparent = True)
    fig2.show()
    fig2.savefig(folder_name + '/Landscape_Kappas.pdf', transparent = True)
    fig2.savefig(folder_name + '/Landscape_Kappas.png', transparent = True)
    fig3.show()
    fig3.savefig(folder_name + '/StableStates_Kappas.pdf', transparent = True)
    fig3.savefig(folder_name + '/StableStates_Kappas.png', transparent = True)