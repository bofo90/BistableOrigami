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

def countStableStates(angles, steps, simulations):
    
    import scipy.cluster.hierarchy as hierarch
    from scipy.spatial.distance import pdist
    
    finalAngles = np.empty((0,np.size(dataAngles,1)))
    for hinge in np.arange(simulations):
        sortAllAngIndex = np.lexsort((angles[steps*(hinge+1)-1,:],angles[steps*hinge,:]))
        finalAngles = np.append(finalAngles, [dataAngles[steps*(hinge+1)-1,sortAllAngIndex]], axis = 0)
                                                         
    Z = hierarch.linkage(finalAngles, 'centroid')
    inverse = hierarch.fcluster(Z, 1, criterion='distance')
    c = hierarch.cophenet(Z, pdist(finalAngles))
    print('this is the cophenet of the hierarchical linkage', c[0])
    
    return inverse

def NiceGraph2D(axes, nameX, nameY, mincoord = [np.NaN, np.NaN], maxcoord = [np.NaN, np.NaN], divisions = [np.NaN, np.NaN],
                buffer = [0.0, 0.0, 0.0]):
    gray = '0.2'
    matl.rcParams.update({'font.size': 11})

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
    return

folder_name = "Results/SingleVertex3/active-set/energy/21-May-2019_temp/kh0.000_kta100.000_ke1.000_kf100.000"

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
flagmask = (exfl !=1) & (exfl !=5) & (exfl !=5)
flagmask = ~flagmask.any(axis= 1)

dataVar['Mask'] = flagmask
dataVar['StableStates'] = countStableStates(dataAngles, IterPerSimulAngles, TotSimul)

dataVar = dataVar.join(dataEnergy[['Hinge Number','TotalEnergy']][IterPerSimulEnergy-2::IterPerSimulEnergy].set_index('Hinge Number'), on = 'HingeNumber')
dataVar.set_index(['Kappa','TargetAngle'], inplace = True)
dataVar.sort_index(inplace=True)

data = dataVar.to_xarray()

fig1 = plt.figure(0,figsize=(cm2inch(24.1), cm2inch(20)))
ax1 = plt.subplot(111)

data['TotalEnergy'].plot(axes = ax1, cmap = matl.cm.nipy_spectral, cbar_kwargs={'ticks': np.linspace(0, data['TotalEnergy'].max(), 5),
    'pad': 0.01})
    
NiceGraph2D(ax1, 'a', 'b')

ax1.xaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g'))
ax1.yaxis.set_major_formatter(matl.ticker.FormatStrFormatter('%.2g'))

#cbar = plt.colorbar(cs1,  ax = ax1, fraction=0.05, pad=0.01, extend = 'max')#format="%d", 
#cbar.set_ticks(np.linspace(0, maxEnergy, 5))
#cbar.set_label('Energy', fontsize = 11, color = '0.2')
#cbar.ax.tick_params(axis='y',colors='0.2')
#cbar.ax.tick_params(axis='x',colors='0.2')
#cbar.outline.set_edgecolor('0.2')